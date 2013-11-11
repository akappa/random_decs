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
#include <chrono>

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/traits.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/operations.hpp>
#include <boost/iostreams/seek.hpp>

namespace {

void print_usage() {
	std::cerr << "Usage:" << std::endl;
	std::cerr << "[-d] infile [-o outfile]" << std::endl;
	std::cerr << "-d\t\tto decompress;" << std::endl;
	std::cerr << "-o outfile\tto select the output filename;" << std::endl;
	std::cerr << "-l level\tcompression level" << std::endl;
}

enum Action {COMPRESS, DECOMPRESS};

struct Arguments {
	std::string input_file;
	std::string output_file;
	unsigned int compress_level;
	Action act;
};

struct File {
	const char *content;
	size_t size;
};

Arguments parse_args(int argc, char **argv) throw (std::runtime_error) {
	Action act = COMPRESS;
	std::string infile;
	std::string outfile;
	unsigned int level = 9;
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
			break;
		case '?':
			throw std::runtime_error("Unknown option: " + optopt);
		}
	}

	if (optind < argc)
		infile = argv[optind];
	else
		throw std::runtime_error("infile not specified");
	if (outfile.size() == 0)
		outfile = infile + ".lz";
	Arguments to_ret = {infile, outfile, level, act};
	return to_ret;
}

File read_file(std::string i_file) {
	std::ifstream file(i_file.c_str(), 
			std::ios_base::in | std::ios_base::binary);
	char *to_put;
	if (!file)
		throw std::runtime_error("Cannot open file");
	// Determine file length
	file.seekg(0, std::ios_base::end);
	size_t length = file.tellg();
	file.seekg(0, std::ios_base::beg);
	// Read it in memory
	to_put = new char[length];
	file.read(to_put, length);
	if (file.fail() || file.bad() 
					|| (file.gcount() != static_cast<std::streamsize>(length)))
		throw std::runtime_error("Cannot read the whole file");
	file.close();
	File to_ret = { to_put , length };
	return to_ret;
}

void write_file(std::string name, File i_file) 
{
	std::ofstream file(name.c_str(), 
			std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
	if (file.good())
		file.write(i_file.content, i_file.size);
	if (!file.good()) {
		file.close();
		throw std::runtime_error("Failure to write to file");
	}
	file.close();
}

template <typename T>
std::streamsize get_size(boost::iostreams::basic_array<T> &array)
{
	std::pair<T*, T*> end = array.output_sequence();
	return end.second - end.first;
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

}

int main(int argc, char **argv)
{
	using namespace boost::iostreams;
	// Step 1: Parse the arguments
	Arguments args;
	try {
		args = parse_args(argc, argv);
	} catch (std::runtime_error &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		print_usage();
		return 1;
	}

	// Step 2: read the input file
	File input = read_file(args.input_file);
	char *mutable_file = const_cast<char*>(input.content);
	basic_array_source<char> input_source(mutable_file, input.size);
	char *output;
	size_t output_file_size = 0;
	count_filter<char> c_filter(&output_file_size);
	// Step 3: compress/decompress
	if (args.act != COMPRESS) {
		const size_t max = 1500000000; // 1.5GM
		output = new char[max];
		basic_array<char> array_sink(output, max);
		filtering_streambuf<boost::iostreams::input> in;
		in.push(zlib_decompressor());
		in.push(input_source);
		filtering_streambuf<boost::iostreams::output> out;
		out.push(c_filter);
		out.push(array_sink);
		auto t1 = std::chrono::high_resolution_clock::now();
		boost::iostreams::copy(in, out);
		auto t2 = std::chrono::high_resolution_clock::now();
		auto spent = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
		std::cout << "Time " << spent.count() << " msecs" << std::endl;
		output_file_size = c_filter.get_count();
	} else {
		const size_t max = 1500000000; // 1.5G
		output = new char[max];
		basic_array<char> array_sink(output, max);
		filtering_streambuf<boost::iostreams::output> out;
		zlib_params compress_params(args.compress_level);
		out.push(zlib_compressor(compress_params));
		out.push(c_filter);
		out.push(array_sink);
		boost::iostreams::copy(input_source, out);
		output_file_size = c_filter.get_count();
		std::cout << "Size " << output_file_size << std::endl;
	}

	File output_f = {output, output_file_size};
	write_file(args.output_file, output_f);
}
