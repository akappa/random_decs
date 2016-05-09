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

void compress(
  std::string infile, std::string outfile, 
  unsigned int quality, unsigned int window, 
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
  bp.lgwin                   = window;
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

  std::cout << "Time " << time_spent << " μs\n"
            << "Size " << sizeof(s_size) + encoded_size << std::endl;

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
    po::positional_options_description pd;
    pd.add("input-file", 1).add("output-file", 1);

    if (op == Operation::COMPRESS) {
      desc.add_options()
          ("quality,q", po::value<unsigned int>()->default_value(11),
           "Compression quality.")
          ("window,w", po::value<unsigned int>()->default_value(22),
           "Window parameter.")
          ("no-dictionary",  "Disables the static dictionary.")
          ("no-relative",    "Disables relative distances.")
          ("no-part",        "Disables block partitioning.")
          ("no-lit-part",    "Disables block partitioning for literals.")
          ("no-len-part",    "Disables block partitioning for insert-and-copy lengths.")
          ("no-dist-part",   "Disables block partitioning for distances.")
          ("no-context",     "Disables context modeling.")
          ("no-context-lit", "Disables context modeling for literals.")
          ("no-context-dst", "Disables context modeling for distances.");

      pd.add("quality", 1);
    } else if (op == Operation::DECOMPRESS) {
      desc.add_options()
        ("print-stats,p", "Print decompression stats on stderr.")
        ("bench-check,b", "Print CSV-separated benchmark check values.");
    }


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
      auto quality        = vm["quality"].as<unsigned int>();
      auto window         = vm["window"].as<unsigned int>();
      auto no_dict        = vm.count("no-dictionary") > 0;
      auto no_rel         = vm.count("no-relative") > 0;
      auto no_part        = vm.count("no-part") > 0;
      auto no_lit_part    = no_part or vm.count("no-lit-part")  > 0;
      auto no_len_part    = no_part or vm.count("no-len-part")  > 0;
      auto no_dist_part   = no_part or vm.count("no-dist-part") > 0;
      auto no_context     = vm.count("no-context") > 0;
      auto no_context_lit = no_context or vm.count("no-context-lit") > 0;
      auto no_context_dst = no_context or vm.count("no-context-dst") > 0;
      compress(infile, outfile, quality, window, no_dict, no_rel, no_lit_part, no_len_part, no_dist_part, no_context_lit, no_context_dst);
    } else {
      auto print_stats    = vm.count("print-stats") > 0;
      auto bench_check    = vm.count("bench-check") > 0;
      decompress(infile, outfile, print_stats, bench_check);
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << "\n"
              << "Command-line options:"  << "\n"
              << desc << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;

}