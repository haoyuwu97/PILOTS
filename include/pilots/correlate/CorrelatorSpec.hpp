#pragma once

#include <cstdint>
#include <cstddef>
#include <string>

#include "pilots/correlate/LagAxis.hpp"

namespace pilots {

// Common correlator configuration parsed from a measure section.
struct CorrelatorSpec {
  std::string type = "exact"; // exact | multitau | fft
  LagAxis axis = LagAxis::Frame;

  double dt = 0.0;           // seconds per timestep
  double timebin_width = 0.0; // seconds (only for axis=timebin)

  std::size_t lag_stride = 1;  // exact + axis=frame only
  std::size_t block_size = 200; // 0 disables SEM

  // timebin post-processing
  std::uint64_t min_pairs_per_bin = 0;
  bool bin_merge = false;

  // multi-tau options
  std::size_t mt_channels = 16;
  std::size_t mt_levels = 10;
};

} // namespace pilots
