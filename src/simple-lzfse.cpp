//-- Internal libraries
#include <io.hpp>
//-- External libraries
#include <boost/program_options.hpp>
//-- Shipped libraries
extern "C" {
  #include <lzfse/lzfse.h>
}
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
  const std::string &in_name, const std::string &out_name
)
{
  // Read file
  size_t insize;
  auto content     = read_file<std::uint8_t>(in_name.c_str(), &insize);
  auto content_ptr = content.get();

  // Allocate output
  auto out_size    = insize * 1.1;
  std::vector<std::uint8_t> output(out_size + sizeof(std::uint32_t));
  auto out_ptr     = output.data();
  auto out_ptr32   = reinterpret_cast<std::uint32_t*>(out_ptr);
  *out_ptr32++     = insize;
  out_ptr          = reinterpret_cast<std::uint8_t*>(out_ptr32);

  // Compress
  std::vector<std::uint8_t> scratch_buffer(lzfse_encode_scratch_size(), 0UL);
  auto t_1         = std::chrono::high_resolution_clock::now();  
  auto comp_size   = lzfse_encode_buffer(out_ptr, out_size, content_ptr, insize, scratch_buffer.data()) + sizeof(*out_ptr32);
  auto t_2         = std::chrono::high_resolution_clock::now();
  auto time_spent  = std::chrono::duration_cast<std::chrono::microseconds>(t_2 - t_1);

  // Write on disk
  std::ofstream out_file(out_name, std::ios_base::binary);
  write_file(out_file, output.data(), comp_size);

  return std::tuple<std::uint64_t, std::uint64_t>(comp_size, time_spent.count());
}

std::uint64_t buffer_decode(
  std::uint8_t *dest_ptr, size_t dest_size, std::uint8_t *src_ptr, size_t src_len
)
{
  std::vector<std::uint8_t> scratch_buffer(lzfse_decode_scratch_size(), 0UL);
  auto t_1 = std::chrono::high_resolution_clock::now();
  lzfse_decode_buffer(
      dest_ptr, dest_size, src_ptr, src_len, scratch_buffer.data()
  );
  auto t_2 = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(t_2 - t_1).count();
}

std::tuple<std::uint64_t, std::uint64_t> decompress(
  const std::string &in_name, const std::string &out_name
)
{
  // Read file
  size_t insize;
  auto content     = read_file<std::uint8_t>(in_name.c_str(), &insize);
  auto content_ptr = content.get();

  // Read uncompressed length
  auto content_ptr_32 = reinterpret_cast<std::uint32_t*>(content_ptr);
  size_t dec_size     = *content_ptr_32++;
  content_ptr         = reinterpret_cast<std::uint8_t*>(content_ptr_32);

  // Allocate output
  std::vector<std::uint8_t> dec_data(dec_size, 0);

  // Decompress
  auto time_spent =  buffer_decode(
      dec_data.data(), dec_size, content_ptr, insize - sizeof(*content_ptr_32)
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
  auto content     = read_file<std::uint8_t>(in_name.c_str(), &insize);
  auto content_ptr = content.get();

 // Read uncompressed length
  auto content_ptr_32 = reinterpret_cast<std::uint32_t*>(content_ptr);
  size_t dec_size     = *content_ptr_32++;
  content_ptr         = reinterpret_cast<std::uint8_t*>(content_ptr_32);

  // Allocate output
  std::vector<std::uint8_t> dec_data(dec_size, 0);

  std::vector<std::uint64_t> dec_times(iterations, 0);

  for (auto &i : dec_times) {
    i = buffer_decode(
      dec_data.data(), dec_size, content_ptr, insize - sizeof(*content_ptr_32)
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

      size_t out_size, time;
      std::tie(out_size, time) = compress(infile, outfile);
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