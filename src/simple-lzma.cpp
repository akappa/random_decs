#include <iostream>
#include <lzma.h>
#include <assert.h>
#include <memory>
#include <string>
#include <stdlib.h>
#include <unistd.h>
#include <stdexcept>
#include <chrono>
#include <io.hpp>

void print_usage() {
	std::cerr << "Usage:" << std::endl;
	std::cerr << "[-l] [-d] infile [-o outfile]" << std::endl;
	std::cerr << "-d\t\tto decompress;" << std::endl;
	std::cerr << "-o outfile\tto select the output filename;" << std::endl;
	std::cerr << "-l level\tcompression level (ignored in decompression)" << std::endl;
}

enum action {COMPRESS, DECOMPRESS};

struct arguments {
	std::string input_file;
	std::string output_file;
	unsigned int compress_level;
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
		outfile = infile + ".rxz";
	arguments to_ret = {infile, outfile, static_cast<unsigned int>(level), act};
	return to_ret;
}

struct lzma_opts {
    lzma_options_lzma options;
    lzma_filter filter_chain[2];
};

void get_preset(uint32_t preset, lzma_opts *opts)
{
    auto options = &(opts->options);
    auto chain = opts->filter_chain;
    lzma_lzma_preset(options, preset);
    chain[0] = {LZMA_FILTER_LZMA2, options};
    chain[1] = {LZMA_VLI_UNKNOWN, 0};
}

int main(int argc, char **argv)
{
    // Parse the command-line
    arguments args;
    try {
        args = parse_args(argc, argv);
    } catch (std::runtime_error e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        print_usage();
        return 1;
    }

    // Read input file and allocate the output buffer
    const char *file_name = args.input_file.c_str();
    std::ifstream file;
    open_file(file, file_name);
    std::streamoff file_len = file_length(file);
    std::unique_ptr<uint8_t[]> data(new uint8_t[file_len]);
    read_data(file, data.get(), file_len);
    size_t out_len;
    std::unique_ptr<uint8_t[]> output;

    typedef uint32_t flen_t;
    if (args.act == COMPRESS) {
        output.reset(new uint8_t[file_len]);
        // Write uncompressed length and preset
        flen_t *r_buf = reinterpret_cast<flen_t*>(output.get());
        *r_buf++ = file_len;
        *r_buf++ = args.compress_level;
        out_len = std::distance(output.get(), reinterpret_cast<uint8_t*>(r_buf));

        // Get the preset
        lzma_opts opts;
        get_preset(args.compress_level, &opts);

        // Encode it
        lzma_ret ret = lzma_raw_buffer_encode(opts.filter_chain, 0, data.get(), file_len, 
                output.get(), &out_len, file_len);
        if (ret != LZMA_OK) {
            std::cerr << "ERROR: Generic error in compression" << std::endl;
            std::cerr << "LZMA error code: " << ret << std::endl;
            return 1;
        }
        std::cout << "Size " << out_len << std::endl;
    } else {
        // Get uncompressed length and preset
        uint32_t preset;
        flen_t *r_buf = reinterpret_cast<flen_t*>(data.get());
        out_len = *r_buf++;
        preset = *r_buf++;

        // Build the preset
        lzma_opts opts;
        get_preset(preset, &opts);

        // Allocate data
        output.reset(new uint8_t[out_len]);
        size_t in_pos = std::distance(data.get(), reinterpret_cast<uint8_t*>(r_buf));

        // Decompress
        size_t dec_pos = 0;
		auto t1 = std::chrono::high_resolution_clock::now();
        lzma_ret ret = lzma_raw_buffer_decode(opts.filter_chain, 0, data.get(), &in_pos, file_len,
                output.get(), &dec_pos, out_len);
		auto t2 = std::chrono::high_resolution_clock::now();
		auto spent = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
        if (ret != LZMA_OK) {
            std::cerr << "ERROR: Generic error in decompression" << std::endl;
            std::cerr << "LZMA error code: " << ret << std::endl;
            return 1;
        }
        assert(dec_pos == out_len);
		std::cout << "Time " << spent.count() << " μs" << std::endl;
    }
    // Write file
    const char *out_file_name = args.output_file.c_str();
    std::ofstream out_file;
    open_file(out_file, out_file_name);
    write_file(out_file, output.get(), out_len);
}
