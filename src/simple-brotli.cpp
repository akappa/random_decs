#include <brotli/encode.h>
#include <brotli/decode.h>

#include <boost/program_options.hpp>

#include <io.hpp>

#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
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

void compress(std::string infile, unsigned int quality, unsigned int window, std::string outfile)
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
  bp.quality = quality;
  bp.lgwin   = window;

  // Write output on file
  size_t encoded_size = theoretical_max;
  auto t_1 = high_resolution_clock::now();
  if (!BrotliCompressBuffer(bp, in_len, data, &encoded_size, out_data)) {
    throw std::logic_error("Error on BrotliCompress");
  }
  auto t_2 = high_resolution_clock::now();
  auto time_spent = duration_cast<microseconds>(t_2 - t_1).count();

  std::cout << "Time " << time_spent << " μs\n"
            << "Size " << sizeof(s_size) + encoded_size << std::endl;

  // Write output on file
  std::ofstream output(outfile, std::ios::binary);
  write_file(output, storage.data(), sizeof(s_size) + encoded_size);
}

void decompress(std::string infile, std::string outfile)
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
  std::cout << "Time " << time_spent << " μs" << std::endl;
  // Write output on file
  std::ofstream output(outfile, std::ios::binary);
  write_file(output, storage.data(), len);
}


enum class Operation { COMPRESS, DECOMPRESS };

Operation parse_tool(const char *name) {
  std::string s_name = name,
              dec    = "decompress",
              enc    = "compress";
  if (s_name == dec) {
    return Operation::DECOMPRESS;
  } else if (s_name == enc) {
    return Operation::COMPRESS;
  }
  throw std::logic_error("Invalid tool " + s_name + " (must be either " + enc + " or " + dec);
}

int main(int argc, char **argv)
{
  namespace po = boost::program_options;
  po::options_description desc;
  po::variables_map vm;

  if (argc < 1 + 1) {
    std::cerr << "ERROR: tool needed (either compress or decompress)" << std::endl;
    return EXIT_FAILURE;
  }

  auto tool = *++argv;
  --argc;

  try {

    auto op = parse_tool(tool);

    desc.add_options()
        ("input-file,i", po::value<std::string>()->required(),
         "Input file.")
        ("output-file,o", po::value<std::string>()->required(),
         "Output file.");

    if (op == Operation::COMPRESS) {
      desc.add_options()
          ("quality,q", po::value<unsigned int>()->default_value(11),
           "Compression quality")
          ("window,w", po::value<unsigned int>()->default_value(22),
           "Window parameter");
    }

    po::positional_options_description pd;
    pd.add("input-file", 1).add("output-file", 1);

    try {
      po::store(po::command_line_parser(argc, argv).options(desc).positional(pd).run(), vm);
      po::notify(vm);
    } catch (boost::program_options::error &e) {
      throw std::runtime_error(e.what());
    }

    // Collect parameters
    auto infile     = vm["input-file"].as<std::string>();
    auto outfile    = vm["output-file"].as<std::string>();

    if (op == Operation::COMPRESS) {
      auto quality = vm["quality"].as<unsigned int>();
      auto window  = vm["window"].as<unsigned int>();
      compress(infile, quality, window, outfile);
    } else {
      decompress(infile, outfile);
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << "\n"
              << "Command-line options:"  << "\n"
              << desc << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;

}