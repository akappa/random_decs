//-- Internal libraries
#include <invalidate_cache.hpp>
#include <io.hpp>
//-- External libraries
#include <boost/program_options.hpp>
#include <lzma.h>
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

using measure_t = std::uint64_t;
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

struct lzma_opts {
  lzma_options_lzma options;
  lzma_filter filter_chain[2];
};

void get_preset(uint32_t preset, lzma_opts *opts)
{
  auto options = &(opts->options);
  auto chain = opts->filter_chain;
  lzma_lzma_preset(options, preset);
  chain[0] = {LZMA_FILTER_LZMA2, options};
  chain[1] = {LZMA_VLI_UNKNOWN, 0};
}

std::tuple<measure_t, measure_t> compress(const std::string &in_name, const std::string &out_name, size_t quality)
{
  // Read input file and allocate the output buffer
  std::ifstream file;
  open_file(file, in_name.c_str());
  std::streamoff file_len = file_length(file);
  std::vector<uint8_t> data(file_len);
  read_data(file, data.data(), file_len);
  size_t out_len;
  std::vector<uint8_t> output(file_len * 1.1 + 200 * (1 << 10)); // 10% max expansion

  // Write uncompressed length and preset
  std::uint32_t *r_buf = reinterpret_cast<std::uint32_t*>(output.data());
  *r_buf++ = file_len;
  *r_buf++ = quality;
  out_len = std::distance(output.data(), reinterpret_cast<uint8_t*>(r_buf));

  // Get the preset
  lzma_opts opts;
  get_preset(quality, &opts);

  // Encode it
  auto t1  = std::chrono::high_resolution_clock::now();
  auto ret = lzma_raw_buffer_encode(
    opts.filter_chain, 0, 
    data.data(),
    file_len, 
    output.data(), &out_len, output.size()
  );
  auto t2 = std::chrono::high_resolution_clock::now();
  auto spent = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
  if (ret != LZMA_OK) {
    std::stringstream ss;
    ss << "ERROR: Generic error in compression; "
       << "LZMA error code: " << ret;
    throw std::logic_error(ss.str());
  }

  // Write file
  std::ofstream out_file;
  open_file(out_file, out_name.c_str());
  write_file(out_file, output.data(), out_len);

  return std::tuple<measure_t, measure_t>(out_len, spent.count());
}

std::tuple<measure_t, measure_t> decompress(const std::string &in_name, const std::string &out_name)
{
  // Read input file
  std::ifstream file;
  open_file(file, in_name.c_str());
  std::streamoff file_len = file_length(file);
  std::vector<uint8_t> data(file_len);
  read_data(file, data.data(), file_len);

  // Get uncompressed length and preset
  std::uint32_t *r_buf = reinterpret_cast<std::uint32_t*>(data.data());
  auto out_len = *r_buf++;
  auto preset = *r_buf++;

  // Allocate output and build the preset
  lzma_opts opts;
  get_preset(preset, &opts);

  // Allocate data
  std::vector<uint8_t> output(out_len);
  size_t in_pos = std::distance(data.data(), reinterpret_cast<uint8_t*>(r_buf));

  // Decompress
  size_t dec_pos = 0;
  auto t1 = std::chrono::high_resolution_clock::now();
  lzma_ret ret = lzma_raw_buffer_decode(opts.filter_chain, 0, data.data(), &in_pos, file_len,
    output.data(), &dec_pos, out_len);
  auto t2 = std::chrono::high_resolution_clock::now();
  auto spent = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
  if (ret != LZMA_OK) {
    std::stringstream ss;
    ss << "ERROR: Generic error in compression; "
       << "LZMA error code: " << ret;
    throw std::logic_error(ss.str());
  }
  assert(dec_pos == out_len);
 
  // Write file
  std::ofstream out_file;
  open_file(out_file, out_name.c_str());
  write_file(out_file, output.data(), out_len);

  return std::tuple<measure_t, measure_t>(out_len, spent.count());
}

std::vector<measure_t> benchmark(const std::string &in_name, size_t tries)
{
  // Read input file
  std::ifstream file;
  open_file(file, in_name.c_str());
  std::streamoff file_len = file_length(file);
  std::vector<uint8_t> data(file_len);
  read_data(file, data.data(), file_len);

  // Get uncompressed length and preset
  std::uint32_t *r_buf = reinterpret_cast<std::uint32_t*>(data.data());
  auto out_len = *r_buf++;
  auto preset = *r_buf++;

  // Allocate output and build the preset
  lzma_opts opts;
  get_preset(preset, &opts);

  // Allocate data
  std::vector<uint8_t> output(out_len);
  size_t in_pos = std::distance(data.data(), reinterpret_cast<uint8_t*>(r_buf));

  // Decompress
  std::vector<measure_t> times(tries, 0);
  for (auto &time : times) {
    auto start_pos = in_pos;
    size_t dec_pos = 0;
    std::fill(output.begin(), output.end(), 0);
    auto t1 = std::chrono::high_resolution_clock::now();
    lzma_ret ret = lzma_raw_buffer_decode(opts.filter_chain, 0, data.data(), &start_pos, file_len,
      output.data(), &dec_pos, out_len);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto spent = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
    if (ret != LZMA_OK) {
      std::stringstream ss;
      ss << "ERROR: Generic error in compression; "
         << "LZMA error code: " << ret;
      throw std::logic_error(ss.str());
    }
    assert(dec_pos == out_len);
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
        auto quality        = vm["quality"].as<unsigned int>();
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