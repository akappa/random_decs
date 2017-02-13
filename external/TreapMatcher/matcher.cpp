#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <tuple>
#include <vector>

#include <matcher/all.hpp>

#include <boost/program_options.hpp>
#include <boost/progress.hpp>

#include <serializer/match_serialize.hpp>

using byte = std::uint8_t;
using std::uint32_t;
using std::vector;
using std::tuple;

vector<uint8_t> read_file(const char *file_name)
{
  std::ifstream t(file_name);
  if (!t) {
    throw std::logic_error(std::string(file_name) + " not found");
  }
  t.seekg(0, std::ios::end);
  size_t size = t.tellg();
  vector<uint8_t> data(size);
  t.seekg(0);
  t.read(reinterpret_cast<char*>(data.data()), size);
  return data;
}

struct comp_params {
  int lgwin;
  int lgblock;
  std::string out_file;
};

template <typename T>
size_t log2(T val)
{
  size_t res = 0;
  while (val > 0) {
    val >>= 1;
    ++res;
  }
  return res;
}

struct match {
  std::uint32_t pos, distance, length;
  match(size_t pos, size_t distance, size_t length) : pos(pos), distance(distance), length(length) { }
};

template <typename os>
os &operator<<(os &stream, const match &m)
{
  stream << "(" << m.pos
         << " <- " << m.pos - m.distance
         << " @ " << m.length
         << ")";
  return stream;
}

class Matcher {
private:
  const comp_params params;
  const size_t block_window;
  size_t input_pos;
  size_t last_processed_pos;
  RingBuffer rb;
  HashToBinaryTree hash;

  tuple<vector<uint32_t>, vector<BackwardMatch>> matches(bool is_last)
  {
    // Allocates output vectors
    auto num_bytes = input_pos - last_processed_pos;
    const size_t max_matches_per_char = 4UL;
    vector<uint32_t> num_matches(num_bytes, 0UL);
    vector<BackwardMatch> matches(num_bytes * max_matches_per_char);

    // End of starting positions for the hasher
    auto position = last_processed_pos;
    const auto store_end =
      (num_bytes >= hash.StoreLookahead())
      ? position + num_bytes - hash.StoreLookahead() + 1
      : position;

    // Begin of input and "input_mask" (?) -- taken from RingBuffer.
    auto input_data = rb.buffer_;
    auto input_mask = rb.mask_;

    // Initialize the hash structure and "stitch" to the previous block (?)
    if (hash.dirty()) {
      hash = HashToBinaryTree(params.lgwin, position, num_bytes, is_last);
    }
    hash.StitchToPreviousBlock(num_bytes, position, input_data, input_mask);

    // Scans every input position.
    size_t cur_match_pos = 0;
    for (auto i = 0; i + hash.HashTypeLength() - 1 < num_bytes; ++i) {

      // Compute current position, maximum distance and maximum length
      const size_t pos = position+ i;
      size_t max_distance = std::min<size_t>(pos, MaxBackwardLimit(params.lgwin));
      size_t max_length = num_bytes - i;
      if (matches.size() < cur_match_pos + HashToBinaryTree::MaxMatches()) {
        matches.resize(cur_match_pos + HashToBinaryTree::MaxMatches());
      }

      // Actually find all matches and store the number of matches
      auto num_found_matches = hash.FindAllMatches(
        input_data, input_mask, pos, max_length, max_distance, matches.data() + cur_match_pos
      );
      num_matches[i] = num_found_matches;

      /* If we have a long match, skip input positions and just store
         that long match; otherwise, go forward as usual. */
      auto cur_match_end = cur_match_pos + num_found_matches;
      if (num_found_matches > 0) {
        const auto match_len = matches[cur_match_end - 1].match_length();
        if (match_len > MAX_ZOPFLI_LEN_QUALITY_11) {
          const auto skip = match_len - 1;
          matches[cur_match_pos++] = matches[cur_match_end - 1];
          num_matches[i] = 1;
          /* Add the tail of the copy to the hasher. */
          hash.StoreHashRange(
            input_data, input_mask, pos + 1, std::min<size_t>(pos + match_len, store_end)
          );
          std::fill(
            std::next(num_matches.begin(), i + 1),
            std::next(num_matches.begin(), i + 1 + skip),
            0UL
          );
          i += skip;
        } else {
          cur_match_pos = cur_match_end;
        }
      }
    }
    return std::make_tuple(num_matches, matches);
  }
public:
  Matcher(comp_params params, size_t block_size)
    : params(params), block_window(1UL << params.lgblock), input_pos(0UL), last_processed_pos(0UL)
  {
    RingBufferInit(&rb);
    RingBufferSetup(params.lgwin, params.lgblock, &rb);
  }

  vector<match> give_data(uint8_t *data, size_t length, bool is_last)
  {
    // Copy input to ring buffer
    if (length > 0) {
      CopyInputToRingBuffer(&rb, length, data);
      input_pos += length;
    }
    if (is_last or (input_pos - last_processed_pos) == block_window) {
      // Perform matching here
      vector<uint32_t> num_matches;
      vector<BackwardMatch> match_vec;
      vector<match> to_ret;
      std::tie(num_matches, match_vec) = matches(is_last);
      auto match_it = match_vec.begin();
      for (auto n : num_matches) {
        for (auto i = 0UL; i < n; ++i) {
          const auto m = *match_it++;
          to_ret.emplace_back(last_processed_pos, m.distance, m.match_length());
        }
        ++last_processed_pos;
      }
      return to_ret;
    }
    return {};
  }
};

class Serialize {
public:
  virtual void serialize(match *begin, match *end) = 0;

  virtual ~Serialize() { }
};

class EmptySerialize : public Serialize {
public:
  EmptySerialize() { }
  void serialize(match *, match *) override { }
};

class FileSerialize : public Serialize {
private:
  std::ofstream file;
  match_serialize::serializer s;
public:
  FileSerialize(const char *file_name, const size_t file_length)
    : file(file_name), s(file, file_length)
  { }

  void serialize(match *begin, match *end) override
  {
    if (begin != end) {
      s.serialize(begin->pos, begin, end);
    }
  }
};


void print_matches(uint8_t *data, size_t file_size, comp_params params)
{
  Matcher match(params, file_size);

  std::unique_ptr<Serialize> s;
  if (params.out_file != "") {
    s = std::make_unique<FileSerialize>(params.out_file.c_str(), file_size);
  } else {
    s = std::make_unique<EmptySerialize>();
  }

  // Process every block of 64KiB
  const size_t buffer_size = 1 << 16; // 64KiB
  boost::progress_display pbar(file_size, std::cerr);
  for (uint8_t * const data_end = data + file_size; data < data_end; data += buffer_size) {
    // Copy data to input buffer
    auto to_copy = std::min<int>(buffer_size, data_end - data);
    auto matches = match.give_data(data, to_copy, data + buffer_size >= data_end);
    auto m_start = matches.data(), m_it = m_start, m_end = m_start + matches.size();
    while (m_it != m_end) {
      if (m_start->pos == m_it->pos) {
        ++m_it;
        continue;
      }
      s->serialize(m_start, m_it);
      m_start = m_it;
    }
    s->serialize(m_start, m_it);
    pbar += to_copy;
  }
}

int main(int argc, char **argv)
{
  using std::string;
  namespace po = boost::program_options;
  po::options_description desc;
  po::variables_map vm;
  try {
    desc.add_options()
        ("input-file,i", po::value<string>()->required(),
         "Input file.")
        ("output-file,o", po::value<string>(),
         "Output file (if needed).")
        ("window-bits,w", po::value<size_t>()->default_value(22),
         "Window size (log_2).")
        ("block-bits,b", po::value<size_t>()->default_value(18),
         "Block size (log_2).");
    po::positional_options_description pd;
    pd.add("input-file", 1).add("output-file", 1);
    try {
      po::store(po::command_line_parser(argc, argv).options(desc).positional(pd).run(), vm);
      po::notify(vm);
    } catch (boost::program_options::error &e) {
      throw std::runtime_error(e.what());
    }


    auto infile     = vm["input-file"].as<string>();
    auto wb         = vm["window-bits"].as<size_t>();
    auto bb         = vm["block-bits"].as<size_t>();
    auto outfile    = std::string("");
    if (vm.count("output-file") > 0) {
      outfile = vm["output-file"].as<string>();
    }

    comp_params params;
    params.lgwin    = wb;
    params.lgblock  = bb;
    params.out_file = outfile;

    auto file_data = read_file(infile.c_str());
    auto t_1 = std::chrono::high_resolution_clock::now();
    print_matches(file_data.data(), file_data.size(), params);
    auto t_2 = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::microseconds>(t_2 - t_1);
    std::cout << "Time\t" << time.count() << "\tμs" << std::endl;

  } catch (std::exception &e) {
    std::cerr << e.what() << "\n"
              << "Command-line options:"  << "\n"
              << desc << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
