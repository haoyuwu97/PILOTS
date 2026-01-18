#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#if PILOTS_HAS_OPENMP
#include <omp.h>
#endif

#include "pilots/app/Runner.hpp"
#include "pilots/config/IniConfig.hpp"
#include "pilots/measures/MeasureRegistry.hpp"

namespace fs = std::filesystem;

namespace {

struct Cli {
  fs::path config;
  int threads = 0; // 0 = use OMP default
  bool list_measures = false;
  bool validate_config = false;
};

void print_usage(const char* argv0) {
  std::cerr
      << "Usage: " << argv0 << " --config <path> [--threads N] [--validate-config]\n"
      << "       " << argv0 << " --list-measures\n"
      << "       " << argv0 << " --version\n";
}

Cli parse_cli(int argc, char** argv) {
  Cli cli;
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    if (a == "--help" || a == "-h") {
      print_usage(argv[0]);
      std::exit(0);
    } else if (a == "--version") {
      std::cout << PILOTS_VERSION_STR << "\n";
      std::exit(0);
    } else if (a == "--config") {
      if (i + 1 >= argc) throw std::runtime_error("--config requires a value");
      cli.config = fs::path(argv[++i]);
    } else if (a == "--threads") {
      if (i + 1 >= argc) throw std::runtime_error("--threads requires a value");
      cli.threads = std::stoi(argv[++i]);
    } else if (a == "--list-measures") {
      cli.list_measures = true;
    } else if (a == "--validate-config") {
      cli.validate_config = true;
    } else {
      throw std::runtime_error("unknown argument: " + a);
    }
  }
  if (!cli.list_measures && cli.config.empty()) {
    throw std::runtime_error("--config is required (or use --list-measures)");
  }
  return cli;
}

} // namespace

int main(int argc, char** argv) {
  try {
    Cli cli = parse_cli(argc, argv);

    if (cli.list_measures) {
      const auto types = pilots::MeasureRegistry::instance().registered_types();
      for (const auto& t : types) {
        std::cout << t << "\n";
      }
      return 0;
    }

#if PILOTS_HAS_OPENMP
    if (cli.threads > 0) {
      omp_set_num_threads(cli.threads);
    }
#endif

    pilots::IniConfig cfg(cli.config);
    pilots::Runner runner(cfg);
    if (cli.validate_config) {
      return runner.validate_config();
    }
    return runner.run();

  } catch (const std::exception& e) {
    std::cerr << "Fatal: " << e.what() << "\n";
    return 1;
  }
}
