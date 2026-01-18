#pragma once

#include <algorithm>
#include <stdexcept>
#include <string>

namespace pilots {

enum class LagAxis {
  Frame,    // lag = \Delta frame index
  Timestep, // lag = \Delta timestep
  TimeBin   // lag = binned \Delta time
};

inline std::string lag_axis_name(LagAxis a) {
  switch (a) {
    case LagAxis::Frame: return "frame";
    case LagAxis::Timestep: return "timestep";
    case LagAxis::TimeBin: return "timebin";
  }
  return "frame";
}

inline LagAxis parse_lag_axis(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return static_cast<char>(::tolower(c)); });
  if (s == "frame" || s == "frames") return LagAxis::Frame;
  if (s == "timestep" || s == "timesteps" || s == "step" || s == "steps") return LagAxis::Timestep;
  if (s == "timebin" || s == "time_bin" || s == "bin" || s == "bins") return LagAxis::TimeBin;
  throw std::runtime_error("invalid lag_axis: '" + s + "' (use frame|timestep|timebin)");
}

} // namespace pilots
