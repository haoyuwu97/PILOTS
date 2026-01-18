#pragma once

#include <cstdint>
#include <vector>

#include "pilots/correlate/LagAxis.hpp"

namespace pilots {

struct CorrelationSeries {
  LagAxis axis = LagAxis::Frame;

  // For axis=frame: lag = lag_frames (integer values stored as double).
  // For axis=timestep: lag = lag_timesteps.
  // For axis=timebin: lag = time_center (seconds).
  std::vector<double> lag;

  // Physical time in seconds for each lag entry.
  // For axis=timebin, time == lag (bin center).
  std::vector<double> time;

  // Correlation/observable value (here: MSD).
  std::vector<double> value;

  // Number of frame-pairs (time origins) contributing to each lag.
  std::vector<std::uint64_t> count_pairs;

  // Standard error of mean estimated from block means (0 if disabled or insufficient blocks).
  std::vector<double> sem;

  // Number of completed blocks used for SEM.
  std::vector<std::uint64_t> n_blocks;

  // Mean timestep difference for each lag (useful for irregular sampling).
  std::vector<double> mean_dtimestep;
};

// Tensor-6 correlation series (xx, yy, zz, xy, xz, yz).
//
// Notes:
// - count_pairs / n_blocks / mean_dtimestep are shared across tensor components.
// - sem_* values are estimated independently per component from block means.
//
// By design, a scalar observable may still be represented by using only the `value_xx`
// channel and leaving other channels as 0. This allows a unified correlator interface.
struct CorrelationSeriesT6 {
  LagAxis axis = LagAxis::Frame;

  std::vector<double> lag;
  std::vector<double> time;

  // Mean values per tensor component.
  std::vector<double> value_xx;
  std::vector<double> value_yy;
  std::vector<double> value_zz;
  std::vector<double> value_xy;
  std::vector<double> value_xz;
  std::vector<double> value_yz;

  // SEM per tensor component.
  std::vector<double> sem_xx;
  std::vector<double> sem_yy;
  std::vector<double> sem_zz;
  std::vector<double> sem_xy;
  std::vector<double> sem_xz;
  std::vector<double> sem_yz;

  std::vector<std::uint64_t> count_pairs;
  std::vector<std::uint64_t> n_blocks;
  std::vector<double> mean_dtimestep;
};

} // namespace pilots
