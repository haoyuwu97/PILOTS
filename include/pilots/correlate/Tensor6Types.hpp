#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace pilots {

// Convention: (xx, yy, zz, xy, xz, yz)
enum Tensor6Index : int {
  T6_XX = 0,
  T6_YY = 1,
  T6_ZZ = 2,
  T6_XY = 3,
  T6_XZ = 4,
  T6_YZ = 5
};

struct Tensor6 {
  std::array<double, 6> v{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double& xx() { return v[T6_XX]; }
  double& yy() { return v[T6_YY]; }
  double& zz() { return v[T6_ZZ]; }
  double& xy() { return v[T6_XY]; }
  double& xz() { return v[T6_XZ]; }
  double& yz() { return v[T6_YZ]; }

  double xx() const { return v[T6_XX]; }
  double yy() const { return v[T6_YY]; }
  double zz() const { return v[T6_ZZ]; }
  double xy() const { return v[T6_XY]; }
  double xz() const { return v[T6_XZ]; }
  double yz() const { return v[T6_YZ]; }
};

// One frame/sample for a selection represented as a 6-channel per-particle field.
// Vector-valued signals (e.g., position, velocity) should be mapped onto (xx,yy,zz)
// and leave (xy,xz,yz) as 0.
struct T6Slot {
  std::int64_t timestep = 0;
  std::vector<double> xx;
  std::vector<double> yy;
  std::vector<double> zz;
  std::vector<double> xy;
  std::vector<double> xz;
  std::vector<double> yz;
};

// Accumulation stats for one lag entry for a tensor observable.
struct LagStatsT6 {
  std::array<double, 6> sum_val{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double sum_dtimestep = 0.0;
  std::uint64_t count = 0;

  // Block statistics (block means over pairs)
  std::array<double, 6> block_cur_sum{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::uint64_t block_cur_count = 0;
  std::array<double, 6> block_mean_sum{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::array<double, 6> block_mean_sumsq{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::uint64_t block_count = 0;
};

} // namespace pilots
