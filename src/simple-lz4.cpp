//-- Internal libraries
#include <io.hpp>
//-- External libraries
#include <boost/program_options.hpp>
//-- Shipped libraries
#include <lz4.h>
#include <lz4hc.h>
//-- Standard libraries
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

enum class Operation { COMPRESS, DECOMPRESS, BENCHMARK };

Operation parse_tool(const char *name) {
  std::string s_name = name,
              dec    = "decompress",
              enc    = "compress",
              bench  = "benchmark";
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
  std::ifstream file;
  open_file(file, in_name.c_str());
  std::streamoff file_len = file_length(file);
  std::vector<char> data(file_len);
  read_data(file, data.data(), file_len);

  const size_t max = file_len * 1.1 + 200 * (1 << 10); // 10% more than input
  std::vector<char> output(max, 0U);
  size_t output_file_size = 0;

  // First, we write the input size on the first 8 bytes
  *reinterpret_cast<std::uint64_t*>(output.data()) = file_len;
  output_file_size = sizeof(std::uint64_t);

  // Then, we write the compress rep.
  auto t1 = std::chrono::high_resolution_clock::now();
  output_file_size += LZ4_compressHC(data.data(), output.data() + output_file_size, file_len);
  auto t2 = std::chrono::high_resolution_clock::now();
  auto spent = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);

  //std::cout << "Writing to " << out_file_name << ", size " << output_file_size << std::endl;
  std::ofstream out_file;
  open_file(out_file, out_name.c_str());
  write_file(out_file, output.data(), output_file_size);

  return std::tuple<measure_t, measure_t>(output_file_size, spent.count());
}

std::tuple<measure_t, measure_t> decompress(const std::string &in_name, const std::string &out_name)
{
  std::ifstream file;
  open_file(file, in_name.c_str());
  std::streamoff file_len = file_length(file);
  std::vector<char> data(file_len);
  read_data(file, data.data(), file_len);

  auto out_size = *reinterpret_cast<std::uint64_t*>(data.data());
  std::vector<char> output(out_size, 0);

  auto t1 = std::chrono::high_resolution_clock::now();
  LZ4_uncompress(data.data() + sizeof(std::uint64_t), output.data(), out_size);
  auto t2 = std::chrono::high_resolution_clock::now();
  auto spent = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
  if (out_size == static_cast<std::uint64_t>(-1)) {
      std::cerr << "ERROR while decompressing" << std::endl;
      exit(1);
  }

  std::ofstream out_file;
  open_file(out_file, out_name.c_str());
  write_file(out_file, output.data(), out_size);

  return std::tuple<measure_t, measure_t>(out_size, spent.count());
}

std::vector<measure_t> benchmark(const std::string &in_name, size_t tries)
{
  std::ifstream file;
  open_file(file, in_name.c_str());
  std::streamoff file_len = file_length(file);
  std::vector<char> data(file_len);
  read_data(file, data.data(), file_len);

  auto out_size = *reinterpret_cast<std::uint64_t*>(data.data());
  std::vector<char> output(out_size, 0);

  std::vector<measure_t> times(tries, 0);

  for (auto &time : times) {
    auto t1 = std::chrono::high_resolution_clock::now();
    LZ4_uncompress(data.data() + sizeof(std::uint64_t), output.data(), out_size);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto spent = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
    if (out_size == static_cast<std::uint64_t>(-1)) {
        std::cerr << "ERROR while decompressing" << std::endl;
        exit(1);
    }
    time = spent.count();    
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
          ("quality,q", po::value<unsigned int>()->default_value(11),
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
