//-- Internal libraries
#include <io.hpp>
//-- External libraries
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/traits.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/operations.hpp>
#include <boost/iostreams/seek.hpp>
#include <boost/program_options.hpp>
//-- Shipped libraries
#include <zopfli/zopfli.h>
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

std::tuple<std::uint64_t, std::uint64_t> compress(
  const std::string &in_name, const std::string &out_name, size_t iterations
)
{
  // Read file
  size_t insize;
  auto content     = read_file<unsigned char>(in_name.c_str(), &insize);
  auto content_ptr = content.get();

  // Set options
  ZopfliOptions opts;
  ZopfliInitOptions(&opts);
  opts.numiterations = iterations;

  // Compress
  unsigned char *out = nullptr;
  size_t comp_size = 0UL;
  auto t_1         = std::chrono::high_resolution_clock::now();  
  ZopfliCompress(
    &opts, ZOPFLI_FORMAT_ZLIB, content_ptr, insize, &out, &comp_size
  );
  auto t_2         = std::chrono::high_resolution_clock::now();
  auto time_spent  = std::chrono::duration_cast<std::chrono::microseconds>(t_2 - t_1);

  // Write on disk
  std::ofstream out_file(out_name, std::ios_base::binary);

  std::uint64_t insize_64 = insize;
  if (!out_file.write(reinterpret_cast<char*>(&insize_64), sizeof(insize_64))) {
    throw std::logic_error("Failed to write on file");
  }

  write_file(out_file, out, comp_size);

  return std::tuple<std::uint64_t, std::uint64_t>(comp_size, time_spent.count());
}

std::uint64_t decompress_buffer(
  const char *in_ptr,
  size_t in_size,
  char *out_ptr,
  size_t out_size
)
{
  using boost::iostreams::basic_array_source;
  using boost::iostreams::basic_array_sink;
  using boost::iostreams::filtering_streambuf;
  using boost::iostreams::zlib_decompressor;
  // Set source
  basic_array_source<char> input_source(in_ptr, in_size);
  filtering_streambuf<boost::iostreams::input> in;
  in.push(zlib_decompressor());
  in.push(input_source);

  // Set sink
  basic_array_sink<char> out(out_ptr, out_size);
  auto t1 = std::chrono::high_resolution_clock::now();
  boost::iostreams::copy(in, out);
  auto t2 = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
}

std::tuple<std::uint64_t, std::uint64_t> decompress(
  const std::string &in_name, const std::string &out_name
)
{
  // Read file
  size_t insize;
  auto content     = read_file(in_name.c_str(), &insize);
  auto content_ptr = content.get();

  // Read uncompressed length
  auto content_ptr_64 = reinterpret_cast<std::uint64_t*>(content_ptr);
  size_t dec_size     = *content_ptr_64++;
  content_ptr         = reinterpret_cast<char*>(content_ptr_64);

  // Allocate output
  std::vector<char> dec_data(dec_size, 0);

  // Decompress
  auto time_spent =  decompress_buffer(
      content_ptr, insize - sizeof(*content_ptr_64), dec_data.data(), dec_size
  );

  // Write on file
  std::ofstream out_file(out_name, std::ios_base::binary);
  write_file(out_file, dec_data.data(), dec_size);

  return std::tuple<std::uint64_t, std::uint64_t>(
    dec_size,
    time_spent
  );
}

std::vector<std::uint64_t> benchmark(
  const std::string &in_name, size_t iterations
)
{
  // Read file
  size_t insize;
  auto content     = read_file<char>(in_name.c_str(), &insize);
  auto content_ptr = content.get();

  // Read uncompressed length
  auto content_ptr_64 = reinterpret_cast<std::uint64_t*>(content_ptr);
  size_t dec_size     = *content_ptr_64++;
  content_ptr         = reinterpret_cast<char*>(content_ptr);

  // Allocate output
  std::vector<char> dec_data(dec_size, 0);

  std::vector<std::uint64_t> dec_times(iterations, 0);

  for (auto &i : dec_times) {
    i = decompress_buffer(
      content_ptr, insize - sizeof(*content_ptr_64), dec_data.data(), dec_size
    );
  }

  return dec_times;
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
             ("iterations,i", po::value<size_t>()->default_value(15),
              "Number of parsing iterations (positive integer).");
        pd.add("iterations", 1);
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
      auto iterations  = vm["iterations"].as<size_t>();

      size_t out_size, time;
      std::tie(out_size, time) = compress(infile, outfile, iterations);
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