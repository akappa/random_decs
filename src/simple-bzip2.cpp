#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <sys/time.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdexcept>
#include "io.hpp"
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/traits.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/operations.hpp>
#include <boost/iostreams/seek.hpp>
#include <chrono>

enum action {COMPRESS, DECOMPRESS};

struct arguments {
  std::string input_file;
  std::string output_file;
  int compress_level;
  action act;
};

struct file {
  const char *content;
  size_t size;
};

arguments parse_args(int argc, char **argv) throw (std::runtime_error) {
  action act = COMPRESS;
  std::string infile;
  std::string outfile;
  int level = 9;
  int c;
  while ((c = getopt(argc, argv, "do:l:")) != -1) {
    switch (c) {
    case 'd':
      act = DECOMPRESS;
      break;
    case 'o':
      outfile = optarg;
      break;
    case 'l':
      level = atoi(optarg);
            if (level < 0 || level > 9)
                throw std::runtime_error("Invalid compression level");
      break;
    case '?':
      throw std::runtime_error("Unknown option: " + std::to_string(optopt));
    }
  }

  if (optind < argc)
    infile = argv[optind];
  else
    throw std::runtime_error("infile not specified");
  if (outfile.size() == 0)
    outfile = infile + ".rb2";
  arguments to_ret = {infile, outfile, static_cast<unsigned int>(level), act};
  return to_ret;
}

void print_usage() {
  std::cerr << "Usage:" << std::endl;
  std::cerr << "[-d] infile [-o outfile]" << std::endl;
  std::cerr << "-d\t\tto decompress;" << std::endl;
  std::cerr << "-o outfile\tto select the output filename;" << std::endl;
  std::cerr << "-l level\tcompression level" << std::endl;
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

int main(int argc, char **argv)
{
  using namespace boost::iostreams;
  // Step 1: Parse the arguments
  arguments args;
  try {
    args = parse_args(argc, argv);
  } catch (std::runtime_error &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    print_usage();
    return 1;
  }

  // Step 2: read the input file
    const char *file_name = args.input_file.c_str();
    std::ifstream file;
    open_file(file, file_name);
    std::streamoff file_len = file_length(file);
    std::unique_ptr<char[]> data(new char[file_len]);
    read_data(file, data.get(), file_len);
  basic_array_source<char> input_source(data.get(), file_len);

    const size_t max = 1500000000; // 1.5G
    std::unique_ptr<char[]> output(new char[max]);
  size_t output_file_size = 0;
  // Step 3: compress/decompress
  if (args.act != COMPRESS) {
        basic_array<char> array_sink(output.get(), max);
        filtering_streambuf<boost::iostreams::input> in;
        in.push(bzip2_decompressor());
        in.push(input_source);
        filtering_streambuf<boost::iostreams::output> out;
        out.push(array_sink);
    auto t1 = std::chrono::high_resolution_clock::now();
        output_file_size = boost::iostreams::copy(in, out);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto spent = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
    std::cout << "Time " << spent.count() << " μs" << std::endl;
  } else {
        basic_array<char> array_sink(output.get(), max);
        count_filter<char> c_filter(&output_file_size);
        filtering_streambuf<boost::iostreams::output> out;
        out.push(bzip2_compressor(args.compress_level));
        out.push(c_filter);
        out.push(array_sink);
        boost::iostreams::copy(input_source, out);
        output_file_size = c_filter.get_count();
        std::cout << "Size " << output_file_size << std::endl;
  }

    const char *out_file_name = args.output_file.c_str();
    std::ofstream out_file;
    open_file(out_file, out_file_name);
    write_file(out_file, output.get(), output_file_size);
}
