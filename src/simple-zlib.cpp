//-- Internal libraries
#include <invalidate_cache.hpp>
#include <io.hpp>
//-- External libraries
#include <boost/program_options.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/traits.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/operations.hpp>
#include <boost/iostreams/seek.hpp>
//-- Standard libraries
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stdint.h>
#include <stdlib.h>
#include <string>
#include <tuple>
#include <unistd.h>

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

template <typename T>
class count_filter : public boost::iostreams::multichar_dual_use_filter {
private:
  size_t *count_ptr;
public:
  count_filter(size_t *ptr) : count_ptr(ptr) {}

  template<typename Source>
  std::streamsize read(Source& src, char* s, std::streamsize n)
  {
    *count_ptr += n;
    return boost::iostreams::read(src, s, n);
  }

  template<typename Sink>
  std::streamsize write(Sink& dest, const char* s, std::streamsize n)
  {
    *count_ptr += n;
    return boost::iostreams::write(dest, s, n);
  }

  template<typename Device>
  void close(Device&, std::ios_base::openmode) 
  { 
  }

  size_t get_count()
  {
    return *count_ptr;
  }
};

using measure_t = std::uint64_t;

std::tuple<measure_t, measure_t> compress(const std::string &in_name, const std::string &out_name, size_t quality)
{
  using boost::iostreams::basic_array_source;
  using boost::iostreams::basic_array;
  using boost::iostreams::filtering_streambuf;
  using boost::iostreams::zlib_compressor;
  using boost::iostreams::zlib_params;

  std::ifstream file;
  open_file(file, in_name.c_str());
  std::streamoff file_len = file_length(file);

  std::vector<char> data(file_len, 0U);
  read_data(file, data.data(), file_len);
  basic_array_source<char> input_source(data.data(), data.data() + file_len);

  const size_t max = (200 * 1<<10) + file_len * 1.2; // 20% more than input
  std::vector<char> output(max, 0U);
  auto out_ptr    = output.data();
  auto out_ptr_64 = reinterpret_cast<std::uint64_t*>(out_ptr);
  *out_ptr_64++   = file_len;
  out_ptr         = reinterpret_cast<char*>(out_ptr_64);

  measure_t size  = 0;
  basic_array<char> array_sink(out_ptr, output.data() + output.size());
  count_filter<char> c_filter(&size);
  filtering_streambuf<boost::iostreams::output> out;
  out.push(zlib_compressor(static_cast<int>(quality)));
  out.push(c_filter);
  out.push(array_sink);
  auto t1        = std::chrono::high_resolution_clock::now();
  boost::iostreams::copy(input_source, out);
  auto t2        = std::chrono::high_resolution_clock::now();
  measure_t time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
  size           = c_filter.get_count() + sizeof(std::uint64_t);

  std::ofstream out_file;
  open_file(out_file, out_name.c_str());
  write_file(out_file, output.data(), size);

  return std::make_tuple(size, time);
}

std::tuple<measure_t, measure_t> decompress(const std::string &in_name, const std::string &out_name)
{
  using boost::iostreams::basic_array_source;
  using boost::iostreams::basic_array;
  using boost::iostreams::filtering_streambuf;
  using boost::iostreams::zlib_decompressor;

  std::ifstream file;
  open_file(file, in_name.c_str());
  std::streamoff file_len = file_length(file);
  std::vector<char> data(file_len);
  read_data(file, data.data(), file_len);
  auto data_ptr = data.data();
  auto data_64  = reinterpret_cast<std::uint64_t*>(data_ptr);
  auto dec_len  = *data_64++;
  data_ptr      = reinterpret_cast<char*>(data_64);

  basic_array_source<char> input_source(data_ptr, file_len - sizeof(std::uint64_t));

  std::vector<char> output(dec_len);
  measure_t size = 0;
  basic_array<char> array_sink(output.data(), dec_len);
  filtering_streambuf<boost::iostreams::input> in;
  in.push(zlib_decompressor());
  in.push(input_source);
  filtering_streambuf<boost::iostreams::output> out;
  out.push(array_sink);
  auto t1 = std::chrono::high_resolution_clock::now();
  size = boost::iostreams::copy(in, out);
  auto t2 = std::chrono::high_resolution_clock::now();
  measure_t time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

  const char *out_file_name = out_name.c_str();
  std::ofstream out_file;
  open_file(out_file, out_file_name);
  write_file(out_file, output.data(), size);

  return std::make_tuple(size, time);
}

std::vector<measure_t> benchmark(const std::string &in_name, size_t tries)
{
  using boost::iostreams::basic_array_source;
  using boost::iostreams::basic_array;
  using boost::iostreams::filtering_streambuf;
  using boost::iostreams::zlib_decompressor;

  std::ifstream file;
  open_file(file, in_name.c_str());
  std::streamoff file_len = file_length(file);
  std::vector<char> data(file_len);
  read_data(file, data.data(), file_len);
  auto data_ptr = data.data();
  auto data_64  = reinterpret_cast<std::uint64_t*>(data_ptr);
  auto dec_len  = *data_64++;
  data_ptr      = reinterpret_cast<char*>(data_64);

  std::vector<measure_t> times(tries, 0);
  for (auto &time : times) {
    basic_array_source<char> input_source(data_ptr, file_len - sizeof(std::uint64_t));

    std::vector<char> output(dec_len);
    basic_array<char> array_sink(output.data(), dec_len);
    filtering_streambuf<boost::iostreams::input> in;
    in.push(zlib_decompressor());
    in.push(input_source);
    filtering_streambuf<boost::iostreams::output> out;
    out.push(array_sink);
    auto t1 = std::chrono::high_resolution_clock::now();
    boost::iostreams::copy(in, out);
    auto t2 = std::chrono::high_resolution_clock::now();
    time = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
    wipe_caches();
  }

  return times;
}

int main(int argc, char **argv)
{
  namespace po = boost::program_options;
  po::options_description desc;
  po::variables_map vm;

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
      desc.add_options()
          ("quality,q", po::value<unsigned int>()->default_value(9),
           "Compression quality.");
      pd.add("quality", 1);
    } else if (op == Operation::BENCHMARK) {
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

    std::vector<measure_t> sizes, times;
    switch (op) {
      case Operation::COMPRESS: {
        auto outfile    = vm["output-file"].as<std::string>();
        auto quality    = vm["quality"].as<unsigned int>();
        measure_t size, time;
        std::tie(size, time) = compress(infile, outfile, quality);
        sizes.push_back(size);
        times.push_back(time);
        break;
      }
      case Operation::DECOMPRESS: {
        auto outfile        = vm["output-file"].as<std::string>();
        measure_t size, time;
        std::tie(size, time) = decompress(infile, outfile);
        sizes.push_back(size);
        times.push_back(time);
        break;
      }
      default: {
        assert(op == Operation::BENCHMARK);
        auto tries = vm["tries"].as<unsigned int>();
        times = benchmark(infile, tries);
      }
    }
    for (auto i : sizes) {
      std::cout << "Size\t" << i << "\tbytes" << std::endl;
    }
    for (auto i : times) {
      std::cout << "Time\t" << i << "\tμs" << std::endl;
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << "\n"
              << "Command-line options:"  << "\n"
              << desc << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}