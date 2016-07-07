//-- Internal libraries
#include <io.hpp>
//-- External libraries
#include <boost/program_options.hpp>
//-- Shipped libraries
#include <lzham.h>
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

lzham_compress_level get_compress_level(size_t level)
{
  switch (level) {
    case 0: return LZHAM_COMP_LEVEL_FASTEST;
    case 1: return LZHAM_COMP_LEVEL_FASTER;
    case 2: return LZHAM_COMP_LEVEL_DEFAULT;
    case 3: return LZHAM_COMP_LEVEL_BETTER;
    case 4: return LZHAM_COMP_LEVEL_UBER;
    default: throw std::logic_error("Invalid compression level");
  }
}

std::tuple<size_t, size_t> compress(
  const std::string &in_name, const std::string &out_name, 
  size_t quality, size_t window, size_t parallelism, 
  bool succinct_c, bool succinct_d
)
{
  // Read file
  size_t insize;
  auto content     = read_file<std::uint8_t>(in_name.c_str(), &insize);
  auto content_ptr = content.get();

  // Set options
  lzham_compress_params params  = {};
  params.m_struct_size          = sizeof(params);
  params.m_dict_size_log2       = window;
  params.m_level                = get_compress_level(quality);
  params.m_max_helper_threads   = parallelism - 1U;
  if (succinct_c) {
    params.m_compress_flags    |= LZHAM_COMP_FLAG_EXTREME_PARSING;
  }
  if (succinct_d) {
    params.m_compress_flags    |= LZHAM_COMP_FLAG_TRADEOFF_DECOMPRESSION_RATE_FOR_COMP_RATIO;
  }

  // Write uncompressed size and compression options
  auto comp_bound  = insize * 1.1 + sizeof(std::uint32_t) * 2;
  std::vector<std::uint8_t> out(comp_bound, 0);
  auto out_ptr     = out.data();
  auto out_sptr    = reinterpret_cast<std::uint32_t*>(out_ptr);
  *out_sptr++      = insize;
  *out_sptr++      = window;
  out_ptr          = reinterpret_cast<std::uint8_t*>(out_sptr);
  
  // Write compressed data
  auto t_1         = std::chrono::high_resolution_clock::now();
  std::uint32_t adler32;
  size_t out_size  = comp_bound - sizeof(std::uint32_t) * 2;
  auto status = lzham_compress_memory(
    &params, out_ptr, &out_size, content_ptr, insize, &adler32
  );
  auto t_2         = std::chrono::high_resolution_clock::now();
  auto time_spent  = std::chrono::duration_cast<std::chrono::microseconds>(t_2 - t_1);
  out_size        += sizeof(*out_sptr) * 2;

  if (status != LZHAM_COMP_STATUS_SUCCESS) {
    throw std::logic_error("Compression failed");
  }


  // Write to file
  std::ofstream out_file(out_name, std::ios_base::binary);
  write_file(out_file, out.data(), out_size);
  return std::tuple<size_t, size_t>(out_size, time_spent.count());
}

std::tuple<size_t, size_t> decompress(const std::string &in_name, const std::string &out_name)
{
  // Read file
  size_t insize;
  auto content      = read_file<std::uint8_t>(in_name.c_str(), &insize);
  auto content_ptr  = content.get();

  // Read decompressed size and comp. window
  auto content_sptr = reinterpret_cast<std::uint32_t*>(content_ptr);
  auto dec_size     = *content_sptr++;
  auto window       = *content_sptr++;
  content_ptr       = reinterpret_cast<std::uint8_t*>(content_sptr);

  // Initialize decompress params
  lzham_decompress_params params = {};
  params.m_struct_size           = sizeof(params);
  params.m_dict_size_log2        = window;
  params.m_decompress_flags      = LZHAM_DECOMP_FLAG_OUTPUT_UNBUFFERED;

  // Allocate output
  std::vector<std::uint8_t> dec_file(dec_size, 0);
  auto dec_ptr      = dec_file.data();

  // Decompress there
  auto t_1         = std::chrono::high_resolution_clock::now();
  size_t d_size_st = dec_size;
  size_t c_size_st = insize - sizeof(*content_sptr) * 2;
  std::uint32_t adler32;
  auto status      = lzham_decompress_memory(
    &params, 
    dec_ptr, &d_size_st, 
    content_ptr, c_size_st, 
    &adler32
  );
  auto t_2         = std::chrono::high_resolution_clock::now();
  auto time_spent  = std::chrono::duration_cast<std::chrono::microseconds>(t_2 - t_1);

  if (status != LZHAM_DECOMP_STATUS_SUCCESS) {
    throw std::logic_error("Decompression failed");
  }

  // Write to file
  std::ofstream out_file(out_name, std::ios_base::binary);
  write_file(out_file, dec_file.data(), dec_size);

  return std::tuple<size_t, size_t>(dec_size, time_spent.count());
}

std::vector<size_t> benchmark(const std::string &in_name, size_t tries)
{
  // Read file
  size_t insize;
  auto content      = read_file<std::uint8_t>(in_name.c_str(), &insize);
  auto content_ptr  = content.get();

  // Read decompressed size and comp. window
  auto content_sptr = reinterpret_cast<std::uint32_t*>(content_ptr);
  auto dec_size     = *content_sptr++;
  auto window       = *content_sptr++;
  content_ptr       = reinterpret_cast<std::uint8_t*>(content_sptr);

  // Initialize decompress params
  lzham_decompress_params params = {};
  params.m_struct_size           = sizeof(params);
  params.m_dict_size_log2        = window;
  params.m_decompress_flags      = LZHAM_DECOMP_FLAG_OUTPUT_UNBUFFERED;

  // Allocate output
  std::vector<std::uint8_t> dec_file(dec_size, 0);
  auto dec_ptr      = dec_file.data();

  size_t d_size_st = dec_size;
  std::uint32_t adler32;
  std::vector<size_t> timings(tries, 0);
  for (auto &i : timings) {
    auto t_1         = std::chrono::high_resolution_clock::now();
    size_t c_size_st = insize - sizeof(*content_sptr) * 2;
    auto status = lzham_decompress_memory(
      &params, 
      dec_ptr, &d_size_st, 
      content_ptr, c_size_st, 
      &adler32
    );
    auto t_2         = std::chrono::high_resolution_clock::now();
    auto time_spent  = std::chrono::duration_cast<std::chrono::microseconds>(t_2 - t_1);
    i                = time_spent.count();
    std::fill(dec_file.begin(), dec_file.end(), 0);
    if (status != LZHAM_DECOMP_STATUS_SUCCESS) {
      throw std::logic_error("Decompression failed");
    }

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
             ("level,l", po::value<size_t>()->default_value(4),
              "Compression level (1-4).")
             ("window,w", po::value<size_t>()->default_value(28),
              "Window size, in log_2(bits): e.g., 28 = 256MiB.")
             ("parallelism,p", po::value<size_t>()->default_value(1),
              "Number of threads used, positive integer.")
             ("succinct-compress,s", po::value<bool>()->default_value(false), "Better compression ratio, slower compression.")
             ("succinct-decompress,S", po::value<bool>()->default_value(false), "Better compression ratio, slower decompression.");
        pd.add("level", 1).add("window", 1);
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
      auto outfile     = vm["output-file"].as<std::string>();
      auto quality     = vm["level"].as<size_t>();
      auto window      = vm["window"].as<size_t>();
      auto parallelism = vm["parallelism"].as<size_t>();
      auto succinct_c  = vm["succinct-compress"].as<bool>();
      auto succinct_d  = vm["succinct-decompress"].as<bool>();

      size_t out_size, time;
      std::tie(out_size, time) = compress(
        infile, outfile, 
        quality, window, parallelism, succinct_c, succinct_d
      );
      std::cout << "Size\t" << out_size << "\tbytes\n"
                << "Time\t" << time     << "\tμs" << std::endl;
    } else if (op == Operation::DECOMPRESS) {
      auto outfile    = vm["output-file"].as<std::string>();
      size_t out_size, time;
      std::tie(out_size, time) = decompress(infile, outfile);
      std::cout << "Size\t" << out_size << "\tbytes\n"
                << "Time\t" << time     << "\tμs" << std::endl;
    } else {
      auto tries  = vm["tries"].as<size_t>();
      auto times  = benchmark(infile, tries);
      for (auto i : times) {
        std::cout << "Time\t" << i << "\tμs" << std::endl;
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