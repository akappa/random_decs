//-- Internal libraries
#include <invalidate_cache.hpp>
#include <io.hpp>
//-- External Shipped Libraries
#include <snappy.h>
#include <snappy-sinksource.h>
//-- External libraries
#include <boost/program_options.hpp>
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
  using std::string;
  string s_name = name,
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

using measure_t = std::uint64_t;

std::tuple<measure_t, measure_t> compress(const std::string &in_name, const std::string &out_name)
{
  // Read input
  size_t in_len;
  auto in_data = read_file<char>(in_name.c_str(), &in_len);
  auto *in_ptr = in_data.get();

  // Allocate output and set sink
  std::vector<char> output(snappy::MaxCompressedLength(in_len));
  snappy::UncheckedByteArraySink sink(output.data());

  // Set source
  snappy::ByteArraySource source(in_ptr, in_len);

  // Compress
  auto t_1      = std::chrono::high_resolution_clock::now();
  auto out_size = snappy::Compress(&source, &sink);
  auto t_2      = std::chrono::high_resolution_clock::now();
  auto spent    = std::chrono::duration_cast<std::chrono::microseconds>(t_2 - t_1);

  // Write output on file
  std::ofstream out_file(out_name, std::ios::binary);
  write_file(out_file, output.data(), out_size);

  return std::tuple<measure_t, measure_t>(out_size, spent.count());
}

std::tuple<measure_t, measure_t> decompress(const std::string &in_name, const std::string &out_name)
{
  // Read input
  size_t in_len;
  auto in_data = read_file<char>(in_name.c_str(), &in_len);
  auto *in_ptr = in_data.get();

  // Allocate output
  size_t out_size;
  snappy::GetUncompressedLength(in_ptr, in_len, &out_size);
  std::vector<char> output(out_size);

  auto t1    = std::chrono::high_resolution_clock::now();
  bool ok    = snappy::RawUncompress(in_ptr, in_len, output.data());
  auto t2    = std::chrono::high_resolution_clock::now();
  auto spent = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);

  if (!ok) {
    throw std::logic_error("Decompression error, exiting");
  }

  // Write output on file
  std::ofstream out_file(out_name, std::ios::binary);
  write_file(out_file, output.data(), out_size);

  return std::tuple<measure_t, measure_t>(out_size, spent.count());
}

std::vector<measure_t> benchmark(const std::string &in_name, size_t tries)
{
  // Read input
  size_t in_len;
  auto in_data = read_file<char>(in_name.c_str(), &in_len);
  auto *in_ptr = in_data.get();

  // Allocate output
  size_t out_size;
  snappy::GetUncompressedLength(in_ptr, in_len, &out_size);
  std::vector<char> output(out_size);

  std::vector<measure_t> times(tries);

  for (auto &time : times) {
    auto t1    = std::chrono::high_resolution_clock::now();
    bool ok    = snappy::RawUncompress(in_ptr, in_len, output.data());
    auto t2    = std::chrono::high_resolution_clock::now();
    auto spent = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);

    if (!ok) {
      throw std::logic_error("Decompression error, exiting");
    }
    time = spent.count();
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
    } else {
      desc.add_options()
        ("tries,d", po::value<size_t>()->default_value(3),
          "Number of times decompression must be performed");
      pd.add("tries",1);
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
        measure_t size, time;
        std::tie(size, time) = compress(infile, outfile);
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
        auto tries = vm["tries"].as<size_t>();
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