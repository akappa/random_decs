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
#include "io.hpp"
#include "lz4.h"
#include "lz4hc.h"

namespace {

enum action {COMPRESS, DECOMPRESS};

struct arguments {
	std::string input_file;
	std::string output_file;
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
	int c;
	while ((c = getopt(argc, argv, "do:l:")) != -1) {
		switch (c) {
		case 'd':
			act = DECOMPRESS;
			break;
		case 'o':
			outfile = optarg;
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
		outfile = infile + ".lz4";
	arguments to_ret = {infile, outfile, act};
	return to_ret;
}

void print_usage() {
	std::cerr << "Usage:" << std::endl;
	std::cerr << "[-d] infile [-o outfile]" << std::endl;
	std::cerr << "-d\t\tto decompress;" << std::endl;
	std::cerr << "-o outfile\tto select the output filename;" << std::endl;
}
}

int main(int argc, char **argv)
{
	// The field is 64-bit long
	typedef std::uint64_t cfile_t;
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

    // Step 3: allocate the output array
    const cfile_t max = 1500000000; // 1.5G
    std::unique_ptr<char[]> output(new char[max]);
	cfile_t output_file_size = 0;

	// Step 4: compress/decompress
	if (args.act != COMPRESS) {
        // First, we retrieve the original data size
        output_file_size = *reinterpret_cast<cfile_t*>(data.get());
		auto t1 = std::chrono::high_resolution_clock::now();
        LZ4_uncompress(data.get() + sizeof(cfile_t), output.get(), output_file_size);
		auto t2 = std::chrono::high_resolution_clock::now();
		auto spent = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        if (output_file_size == static_cast<cfile_t>(-1)) {
            std::cerr << "ERROR while decompressing" << std::endl;
            exit(1);
        }
		std::cout << "Time " << spent.count() << " msecs" << std::endl;
	} else {
        // First, we write the input size on the first 8 bytes
        *reinterpret_cast<cfile_t*>(output.get()) = file_len;
        output_file_size = sizeof(cfile_t);
        // Then, we write the compress rep.
        output_file_size += LZ4_compressHC(data.get(), output.get() + output_file_size, file_len);
        std::cout << "Size " << output_file_size << std::endl;
	}

    const char *out_file_name = args.output_file.c_str();
    //std::cout << "Writing to " << out_file_name << ", size " << output_file_size << std::endl;
    std::ofstream out_file;
    open_file(out_file, out_file_name);
    write_file(out_file, output.get(), output_file_size);
}
