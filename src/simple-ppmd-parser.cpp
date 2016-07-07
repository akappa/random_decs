#include <boost/program_options.hpp>
#include <iostream>

enum class Operation { COMPRESS, DECOMPRESS, BENCHMARK };

Operation parse_tool(const char *name) {
  std::string s_name = name,
              dec    = "decompress",
              enc    = "compress",
              bench = "benchmark";
  if (s_name == dec) {
    return Operation::DECOMPRESS;
  } else if (s_name == enc) {
    return Operation::COMPRESS;
  } else if (s_name == bench) {
    return Operation::BENCHMARK;
  }
  throw std::logic_error("Invalid tool " + s_name + " (must be either " + enc + ", " + dec + " or " + bench);
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
    }

    if (op == Operation::COMPRESS) {
      desc.add_options()
          ("model-order,M", po::value<unsigned int>()->default_value(8),
           "PPMd model order.")
          ("working-memory,m", po::value<unsigned int>()->default_value(128),
           "PPMd working memory, in megabytes. Maximum value: 2048.");
      pd.add("working-memory", 1);
    }else {
      desc.add_options()
        ("tries,t", po::value<unsigned int>()->default_value(3),
          "Number of decompressions");
      pd.add("tries", 1);
    }

    try {
      po::store(po::command_line_parser(argc, argv).options(desc).positional(pd).run(), vm);
      po::notify(vm);
    } catch (boost::program_options::error &e) {
      throw std::runtime_error(e.what());
    }

    // Collect parameters
    auto infile     = vm["input-file"].as<std::string>();
    switch (op) {
      case Operation::COMPRESS: {
        auto outfile     = vm["output-file"].as<std::string>();
        auto model_order = vm["model-order"].as<unsigned int>();
        auto wm          = vm["working-memory"].as<unsigned int>();
        std::cout << "C"         << "\n"
                  << infile      << "\n"
                  << outfile     << "\n"
                  << model_order << "\n"
                  << wm          << std::endl;
        break;
      }
      case Operation::DECOMPRESS: {
        auto outfile        = vm["output-file"].as<std::string>();
        std::cout << "D"         << "\n"
                  << infile      << "\n"
                  << outfile     << std::endl;
        break;
      }
      default: {
        assert(op == Operation::BENCHMARK);
        auto tries = vm["tries"].as<unsigned int>();
        std::cout << "B"         << "\n"
                  << infile      << "\n"
                  << tries       << std::endl;
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