#include <brotli/encode.h>
#include <brotli/decode.h>

#include <invalidate_cache.hpp>

#include <boost/program_options.hpp>

#include <io.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <sstream>
#include <string>

template <typename T>
T read(std::uint8_t *data)
{
  return *reinterpret_cast<T*>(data);
}

template <typename T>
std::uint8_t *write(std::uint8_t *data, T val)
{
  *reinterpret_cast<T*>(data) = val;
  return data + sizeof(val);
}

using byte    = std::uint8_t;
using s_size  = std::uint64_t;

class matcher_option {
public:
  virtual void specific_options(boost::program_options::options_description &options) = 0;
  virtual std::unique_ptr<brotli::ext::matcher> instantiate(boost::program_options::variables_map &vm) = 0;
  virtual std::string option_name() const= 0;
  virtual ~matcher_option() { }
};

struct treap_option : public matcher_option{

  treap_option() {}

  void specific_options(boost::program_options::options_description &options) override
  {
    options.add_options()(
      "window,w", boost::program_options::value<unsigned int>()->default_value(22),
      "Window parameter."
    )(
      "enable-static-dictionary", "Enable Brotli's static dictionary."
    );
  }
  std::unique_ptr<brotli::ext::matcher> instantiate(boost::program_options::variables_map &vm) override
  {
    auto window = vm["window"].as<unsigned int>();
    return std::make_unique<ext::treap>(window);
  }

  std::string option_name() const override
  {
    return "treap";
  }
};

void compress(
  std::string infile, std::string outfile, 
  unsigned int quality, ext::matcher *matcher,
  bool no_dict, bool no_rel,
  bool no_lit_part, bool no_len_part, bool no_dist_part, 
  bool no_context_literal, bool no_context_distance
)
{

  using namespace std::chrono;

  // Read input
  std::ifstream input(infile, std::ios::binary);
  size_t in_len;
  auto in_data = read_file<byte>(infile.c_str(), &in_len);
  auto *data = in_data.get();

  // Figure out maximum theoretical comp. size
  size_t theoretical_max = 3 * in_len + 500;
  
  // Allocate output
  std::vector<byte> storage(theoretical_max, byte{});
  byte *out_data = storage.data();

  // Write uncompressed size
  out_data = write<s_size>(out_data, in_len);

  // Compress (in-memory)
  brotli::BrotliParams bp;
  bp.quality                 = quality;
  bp.matcher                 = matcher;
  bp.enable_dictionary       = !no_dict;
  bp.enable_relative         = !no_rel;
  bp.enable_lit_part         = !no_lit_part;
  bp.enable_len_part         = !no_len_part;
  bp.enable_dist_part        = !no_dist_part;
  bp.enable_context_literal  = !no_context_literal;
  bp.enable_context_distance = !no_context_distance;

  // Write output on file
  size_t encoded_size = theoretical_max;
  auto t_1 = high_resolution_clock::now();
  if (!BrotliCompressBuffer(bp, in_len, data, &encoded_size, out_data)) {
    throw std::logic_error("Error on BrotliCompress");
  }
  auto t_2 = high_resolution_clock::now();
  auto time_spent = duration_cast<microseconds>(t_2 - t_1).count();

  std::cout << "Time\t" << time_spent << "\tμs\n"
            << "Size\t" << sizeof(s_size) + encoded_size << "\tbytes" << std::endl;

  // Write output on file
  std::ofstream output(outfile, std::ios::binary);
  write_file(output, storage.data(), sizeof(s_size) + encoded_size);
}

void decompress(std::string infile, std::string outfile, bool print_stats, bool bench_check)
{
  using namespace std::chrono;
  // Read input
  std::ifstream input(infile, std::ios::binary);
  size_t in_len;
  auto in_data = read_file<byte>(infile.c_str(), &in_len);
  auto *data = in_data.get();

  // Figure out uncompressed size
  auto len   = read<s_size>(data);
  data      += sizeof(s_size);
  in_len    -= sizeof(s_size);

  // Alocate output
  std::vector<byte> storage(len, byte{});

  // Decompress (in-memory)
  size_t dec_len = len;
  BrotliEnableStatsPrint = print_stats;
  BrotliEnableBenchPrint = bench_check;
  auto t_1 = high_resolution_clock::now();
  auto res = BrotliDecompressBuffer(in_len, data, &dec_len, storage.data());
  auto t_2 = high_resolution_clock::now();

  if (res != BROTLI_RESULT_SUCCESS) {
   std::stringstream ss;
   ss << "Decompression returned error code " << res;
   throw std::logic_error(ss.str());
  }
  
  if (len != dec_len) {
   std::cerr << "WARNING: len != dec_len (" << len << " != " << dec_len << ")" << std::endl;
  }

  auto time_spent = duration_cast<microseconds>(t_2 - t_1).count();
  std::cout << "Time\t"  << time_spent << "\tμs\n"
            << "Size\t"  << len        << "\tbytes" << std::endl;
  // Write output on file
  std::ofstream output(outfile, std::ios::binary);
  write_file(output, storage.data(), len);
}

void benchmark(std::string infile, size_t tries)
{
  using namespace std::chrono;
  // Read input
  std::ifstream input(infile, std::ios::binary);
  size_t in_len;
  auto in_data = read_file<byte>(infile.c_str(), &in_len);
  auto *data = in_data.get();

  // Figure out uncompressed size
  auto len   = read<s_size>(data);
  data      += sizeof(s_size);
  in_len    -= sizeof(s_size);

  // Alocate output
  std::vector<byte> storage(len, byte{});

  // Decompress (in-memory)
  std::vector<std::uint64_t> times(tries, 0UL);
  size_t dec_len = len;
  for (auto &i : times) {
    auto t_1 = high_resolution_clock::now();
    auto res = BrotliDecompressBuffer(in_len, data, &dec_len, storage.data());
    auto t_2 = high_resolution_clock::now();
    wipe_caches();
    if (res != BROTLI_RESULT_SUCCESS) {
     std::stringstream ss;
     ss << "Decompression returned error code " << res;
     throw std::logic_error(ss.str());
    }
    if (len != dec_len) {
      std::stringstream ss;

      ss << "WARNING: len != dec_len (" << len << " != " << dec_len << ")";
      throw std::logic_error(ss.str());
    }
    i          = std::chrono::duration_cast<std::chrono::microseconds>(t_2 - t_1).count();    
  }

  for (auto i : times) {
    std::cout << "Time\t" << i << "\tμs" << std::endl;
  }
}

enum class Operation { COMPRESS, DECOMPRESS, BENCHMARK };

Operation parse_tool(const char *name) {
  std::string s_name = name,
              dec    = "decompress",
              enc    = "compress",
              bench = "benchmark";
  if (s_name == dec) {
    return Operation::DECOMPRESS;
  } else if (s_name == enc) {
    return Operation::COMPRESS;
  } else if (s_name == bench) {
    return Operation::BENCHMARK;
  }
  throw std::logic_error("Invalid tool " + s_name + " (must be either " + enc + ", " + dec + " or " + bench);
}

int main(int argc, char **argv)
{
  using std::string;
  namespace po = boost::program_options;
  po::options_description desc;
  po::variables_map vm;
  treap_option treap;
  std::vector<matcher_option*> options {{ &treap }};
  auto get_option = [&options] (const std::string &option_name) -> matcher_option* {
    auto opt = std::find_if(
      options.begin(), options.end(),
      [&option_name] (const matcher_option *i) { return i->option_name() ==  option_name; }
    );
    if (opt == options.end()) {
      throw std::logic_error("invalid parser");
    }
    return *opt;
  };


  if (argc < 1 + 1) {
    std::cerr << "ERROR: tool needed (either compress, decompress or benchmark)" << std::endl;
    return EXIT_FAILURE;
  }

  auto tool = *++argv;
  --argc;

  try {

    auto op = parse_tool(tool);

    desc.add_options()
        ("input-file,i", po::value<std::string>()->required(),
         "Input file.");
    po::positional_options_description pd;
    pd.add("input-file", 1);

    if (op != Operation::BENCHMARK) {
        desc.add_options()
          ("output-file,o", po::value<std::string>()->required(),
           "Output file.");
        pd.add("output-file", 1);
    }

    if (op == Operation::COMPRESS) {
      std::stringstream ss;
      ss << "Matchers. Alternatives: ";
      for (auto &i : options) {
        ss << i->option_name() << " ";
      }
      auto matcher_options = ss.str();

      desc.add_options()
          ("quality,q", po::value<unsigned int>()->default_value(11),
           "Compression quality.")
          ("matcher,m",   po::value<string>()->default_value(options.front()->option_name()), matcher_options.c_str())
          ("dictionary",  po::value<bool>()->default_value(true), "Enables the static dictionary.")
          ("relative",    po::value<bool>()->default_value(true), "Enables relative distances.")
          ("part",        po::value<bool>()->default_value(true), "Enables block partitioning.")
          ("lit-part",    po::value<bool>()->default_value(true), "Enables block partitioning for literals.")
          ("len-part",    po::value<bool>()->default_value(true), "Enables block partitioning for insert-and-copy lengths.")
          ("dist-part",   po::value<bool>()->default_value(true), "Enables block partitioning for distances.")
          ("context",     po::value<bool>()->default_value(true), "Enables context modeling.")
          ("context-lit", po::value<bool>()->default_value(true), "Enables context modeling for literals.")
          ("context-dst", po::value<bool>()->default_value(true), "Enables context modeling for distances.");
      pd.add("quality", 1);
    } else if (op == Operation::DECOMPRESS) {
      desc.add_options()
        ("print-stats,p", po::value<bool>()->default_value(false), "Print decompression stats on stderr.")
        ("bench-check,b", po::value<bool>()->default_value(false), "Print CSV-separated benchmark check values.");
    } else {
      desc.add_options()
        ("tries,t", po::value<unsigned int>()->default_value(3),
          "Number of decompressions");
      pd.add("tries", 1);
    }

    try {
      po::store(po::command_line_parser(argc, argv).options(desc).positional(pd).run(), vm);
      po::notify(vm);
    } catch (boost::program_options::error &e) {
      throw std::runtime_error(e.what());
    }

    // Collect parameters
    auto infile     = vm["input-file"].as<std::string>();

    switch (op) {
      case Operation::COMPRESS: {
        auto outfile        = vm["output-file"].as<string>();
        auto quality        = vm["quality"].as<unsigned int>();
        auto no_dict        = not vm["dictionary"].as<bool>();
        auto no_rel         = not vm["relative"].as<bool>();
        auto no_part        = not vm["part"].as<bool>();
        auto no_lit_part    = no_part or (not vm["lit-part"].as<bool>());
        auto no_len_part    = no_part or (not vm["len-part"].as<bool>());
        auto no_dist_part   = no_part or (not vm["dist-part"].as<bool>());
        auto no_context     = not vm["context"].as<bool>();
        auto no_context_lit = no_context or (not vm["context-lit"].as<bool>());
        auto no_context_dst = no_context or (not vm["context-dst"].as<bool>());

        auto matcher = get_option(vm["matcher"].as<string>())->instantiate(vm);
        
        compress(infile, outfile, quality, matcher.get(), no_dict, no_rel, no_lit_part, no_len_part, no_dist_part, no_context_lit, no_context_dst);
        break;
      }
      case Operation::DECOMPRESS: {
        auto outfile        = vm["output-file"].as<std::string>();
        auto print_stats    = vm["print-stats"].as<bool>();
        auto bench_check    = vm["bench-check"].as<bool>();
        decompress(infile, outfile, print_stats, bench_check);
        break;
      }
      default: {
        assert(op == Operation::BENCHMARK);
        auto tries = vm["tries"].as<unsigned int>();
        benchmark(infile, tries);
      }
    }
  } catch (std::exception &e) {
    std::cerr << "ERROR: "
              << e.what() << ".\n\n"
              << "Command-line options:"  << "\n"
              << desc << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;

}
