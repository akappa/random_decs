#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <sys/time.h>
#include <stdint.h>
#include <stdexcept>
#include <chrono>

#include <snappy.h>
#include <snappy-sinksource.h>

namespace {

void print_usage() {
	std::cerr << "Usage:" << std::endl;
	std::cerr << "[-d] infile [-o outfile]" << std::endl;
	std::cerr << "-d\t\tto decompress;" << std::endl;
	std::cerr << "-o outfile\tto select the output filename;" << std::endl;
}

enum Action {COMPRESS, DECOMPRESS};

struct Arguments {
	std::string input_file;
	std::string output_file;
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
	int c;
	while ((c = getopt(argc, argv, "do:")) != -1) {
		switch (c) {
		case 'd':
			act = DECOMPRESS;
			break;
		case 'o':
			outfile = optarg;
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
		outfile = infile + ".snp";
	Arguments to_ret = {infile, outfile, act};
	return to_ret;
}

File read_file(std::string i_file) {
	std::ifstream file(i_file.c_str(), std::ios_base::in | std::ios_base::binary);
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
	if (file.fail() || file.bad() || (file.gcount() != static_cast<std::streamsize>(length)))
		throw std::runtime_error("Cannot read the whole file");
	file.close();
	File to_ret = { to_put , length };
	return to_ret;
}

void write_file(std::string name, File i_file) {
	std::ofstream file(name.c_str(), std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
	if (file.good())
		file.write(i_file.content, i_file.size);
	if (!file.good()) {
		file.close();
		throw std::runtime_error("Failure to write to file");
	}
	file.close();
}

}

int main(int argc, char **argv) {
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

	// Step 3: compress/decompress
	char *str_output;
	size_t out_size;
	if (args.act == COMPRESS) {
		str_output = new char[snappy::MaxCompressedLength(input.size)];
		snappy::UncheckedByteArraySink sink(str_output);
		snappy::ByteArraySource source(input.content, input.size);
		out_size = snappy::Compress(&source, &sink);
		std::cout << "Size " << out_size << std::endl;
	} else {
		snappy::GetUncompressedLength(input.content, input.size, &out_size);
		str_output = new char[out_size];
		auto t1 = std::chrono::high_resolution_clock::now();
		bool ok = snappy::RawUncompress(input.content, input.size, str_output);
		auto t2 = std::chrono::high_resolution_clock::now();
		auto spent = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);

		if (!ok) {
			std::cerr << "Decompression error, exiting" << std::endl;
			return 1;
		}
		std::cout << "Time " << spent.count() << " msecs" << std::endl;
	}

	// Step 4: write the output file
	File output_f = { str_output, out_size };
	write_file(args.output_file, output_f);
	return 0;
}
