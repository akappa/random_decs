//-- Internal libraries
#include <io.hpp>
//-- External libraries
#include <boost/program_options.hpp>
//-- Shipped external libraries
#include <zstd/zstd.h>
//-- Standard libraries
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <tuple>

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

std::tuple<size_t, size_t> compress(const std::string &in_name, const std::string &out_name, size_t quality)
{
  // Read file
  size_t insize;
  auto content     = read_file(in_name.c_str(), &insize);
  auto content_ptr = content.get();
  auto comp_bound  = ZSTD_compressBound(insize) + sizeof(std::uint64_t);
  std::vector<char> out(comp_bound, 0);
  auto out_ptr     = out.data();
  // Write uncompressed size
  auto out_sptr    = reinterpret_cast<std::uint64_t*>(out_ptr);
  *out_sptr++      = insize;
  out_ptr          = reinterpret_cast<char*>(out_sptr);
  // Write compressed size
  auto t_1         = std::chrono::high_resolution_clock::now();
  auto out_size    = ZSTD_compress(out_ptr, comp_bound, content_ptr, insize, quality) + sizeof(*out_sptr);
  auto t_2         = std::chrono::high_resolution_clock::now();
  auto time_spent  = std::chrono::duration_cast<std::chrono::microseconds>(t_2 - t_1);
  // Write to file
  std::ofstream out_file(out_name, std::ios_base::binary);
  write_file(out_file, out.data(), out_size);
  return std::tuple<size_t, size_t>(out_size, time_spent.count());
}

std::tuple<size_t, size_t> decompress(const std::string &in_name, const std::string &out_name)
{
  // Read file
  size_t insize;
  auto content      = read_file(in_name.c_str(), &insize);
  auto content_ptr  = content.get();
  // Read decompressed size
  auto content_sptr = reinterpret_cast<std::uint64_t*>(content_ptr);
  auto dec_size     = *content_sptr++;
  content_ptr       = reinterpret_cast<char*>(content_sptr);
  // Allocate output
  std::vector<char> dec_file(dec_size, 0);
  auto dec_ptr      = dec_file.data();
  // Decompress there
  auto t_1         = std::chrono::high_resolution_clock::now();
  ZSTD_decompress(dec_ptr, dec_size, content_ptr, insize - sizeof(content_sptr));
  auto t_2         = std::chrono::high_resolution_clock::now();
  auto time_spent  = std::chrono::duration_cast<std::chrono::microseconds>(t_2 - t_1);
  // Write to file
  std::ofstream out_file(out_name, std::ios_base::binary);
  write_file(out_file, dec_file.data(), dec_size);

  return std::tuple<size_t, size_t>(dec_size, time_spent.count());
}

std::vector<size_t> benchmark(const std::string &in_name, size_t tries)
{
  // Read file
  size_t insize;
  auto content      = read_file(in_name.c_str(), &insize);
  auto content_ptr  = content.get();
  // Read decompressed size
  auto content_sptr = reinterpret_cast<std::uint64_t*>(content_ptr);
  auto dec_size     = *content_sptr++;
  content_ptr       = reinterpret_cast<char*>(content_sptr);
  // Allocate output
  std::vector<char> dec_file(dec_size, 0);
  auto dec_ptr      = dec_file.data();

  std::vector<size_t> timings(tries, 0);
  for (auto &i : timings) {
    auto t_1         = std::chrono::high_resolution_clock::now();
    ZSTD_decompress(dec_ptr, dec_size, content_ptr, insize - sizeof(content_sptr));
    auto t_2         = std::chrono::high_resolution_clock::now();
    auto time_spent  = std::chrono::duration_cast<std::chrono::microseconds>(t_2 - t_1);
    i                = time_spent.count();
    std::fill(dec_file.begin(), dec_file.end(), 0);
  }
  return timings;
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
      if (op == Operation::COMPRESS) {
        desc.add_options()
             ("level", po::value<size_t>()->default_value(22),
              "Compression level (1-22, default: 22).");
        pd.add("level", 1);
      }
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

    if (op == Operation::COMPRESS) {
      auto outfile    = vm["output-file"].as<std::string>();
      auto quality    = vm["level"].as<size_t>();
      size_t out_size, time;
      std::tie(out_size, time) = compress(infile, outfile, quality);
      std::cout << "Size\t" << out_size << "\n"
                << "Time\t" << time     << std::endl;
    } else if (op == Operation::DECOMPRESS) {
      auto outfile    = vm["output-file"].as<std::string>();
      size_t out_size, time;
      std::tie(out_size, time) = decompress(infile, outfile);
      std::cout << "Size\t" << out_size << "\n"
                << "Time\t" << time     << std::endl;
    } else {
      auto tries  = vm["tries"].as<size_t>();
      auto times  = benchmark(infile, tries);
      for (auto i : times) {
        std::cout << "Time\t" << i << std::endl;
      }
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << "\n"
              << "Command-line options:"  << "\n"
              << desc << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}