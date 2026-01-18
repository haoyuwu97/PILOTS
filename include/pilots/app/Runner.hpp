#pragma once

#include <filesystem>

#include "pilots/config/IniConfig.hpp"

namespace pilots {

// Runner: ...
//
// P1-B core goal:
//   main() only handles CLI + config, then calls Runner(cfg).run().
//   Runner owns the pipeline: on_start -> on_frame -> flush/checkpoint -> finalize.
class Runner {
public:
  explicit Runner(const IniConfig& cfg);

  // Execute the run. Returns 0 on success.
  int run();

  // Validate config + input dependencies and exit without processing the full
  // trajectory (CLI: --validate-config). This performs the same fail-fast
  // checks as run() up to and including measure instantiation and centralized
  // dependency validation, but skips the frame-processing loop and avoids
  // writing results.
  int validate_config();

private:
  const IniConfig& cfg_;

  int run_impl_(bool validate_only);
};

} // namespace pilots
