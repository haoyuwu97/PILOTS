#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "pilots/core/Frame.hpp"
#include "pilots/select/Selection.hpp"

namespace pilots {

class MSD {
public:
  MSD() = default;

  // Choose which Cartesian components contribute to the MSD.
  // Supported strings in config: "xyz" (default), "x", "y", "z", "xy", "xz", "yz".
  void set_component_mask(int mask) {
    if (mask <= 0 || mask > 7) {
      throw std::runtime_error("MSD: invalid component mask (must be 1..7)");
    }
    comp_mask_ = mask;
  }

  void set_origin(const Frame& frame) {
    const std::size_t n = frame.natoms;
    if (n == 0) throw std::runtime_error("MSD: cannot set origin for empty frame");
    x0_.assign(frame.xu.begin(), frame.xu.end());
    y0_.assign(frame.yu.begin(), frame.yu.end());
    z0_.assign(frame.zu.begin(), frame.zu.end());
    origin_timestep_ = frame.timestep;
    initialized_ = true;
  }

  bool initialized() const { return initialized_; }
  std::int64_t origin_timestep() const { return origin_timestep_; }

  // Compute MSD over `sel` indices.
  // If remove_drift is true, subtract the mean displacement vector computed over `drift_sel`:
  //   MSD = <|dr - v|^2>_sel, where v = <dr>_drift
  double compute(const Frame& frame,
                 const Selection& sel,
                 bool remove_drift,
                 const Selection& drift_sel) const {
    if (!initialized_) throw std::runtime_error("MSD: origin not set");
    if (frame.natoms != x0_.size()) throw std::runtime_error("MSD: natoms mismatch");
    if (sel.idx.empty()) throw std::runtime_error("MSD: selection '" + sel.name + "' is empty");

    const std::size_t nsel = sel.idx.size();

    const bool use_x = (comp_mask_ & 1) != 0;
    const bool use_y = (comp_mask_ & 2) != 0;
    const bool use_z = (comp_mask_ & 4) != 0;

    double sum_dr2 = 0.0;
    double sum_dx_sel = 0.0, sum_dy_sel = 0.0, sum_dz_sel = 0.0;

#if PILOTS_HAS_OPENMP
#pragma omp parallel for reduction(+:sum_dr2,sum_dx_sel,sum_dy_sel,sum_dz_sel)
#endif
    for (std::size_t k = 0; k < nsel; ++k) {
      const std::size_t i = sel.idx[k];
      const double dx = frame.xu[i] - x0_[i];
      const double dy = frame.yu[i] - y0_[i];
      const double dz = frame.zu[i] - z0_[i];
      if (use_x) sum_dx_sel += dx;
      if (use_y) sum_dy_sel += dy;
      if (use_z) sum_dz_sel += dz;
      double dr2 = 0.0;
      if (use_x) dr2 += dx*dx;
      if (use_y) dr2 += dy*dy;
      if (use_z) dr2 += dz*dz;
      sum_dr2 += dr2;
    }

    const double inv_sel = 1.0 / static_cast<double>(nsel);
    const double mean_dr2_sel = sum_dr2 * inv_sel;

    if (!remove_drift) {
      return mean_dr2_sel;
    }

    if (drift_sel.idx.empty()) {
      throw std::runtime_error("MSD: drift selection is empty");
    }

    // compute mean drift vector v over drift_sel
    const std::size_t nd = drift_sel.idx.size();
    double sum_dx_d = 0.0, sum_dy_d = 0.0, sum_dz_d = 0.0;

#if PILOTS_HAS_OPENMP
#pragma omp parallel for reduction(+:sum_dx_d,sum_dy_d,sum_dz_d)
#endif
    for (std::size_t k = 0; k < nd; ++k) {
      const std::size_t i = drift_sel.idx[k];
      const double dx = frame.xu[i] - x0_[i];
      const double dy = frame.yu[i] - y0_[i];
      const double dz = frame.zu[i] - z0_[i];
      if (use_x) sum_dx_d += dx;
      if (use_y) sum_dy_d += dy;
      if (use_z) sum_dz_d += dz;
    }

    const double inv_d = 1.0 / static_cast<double>(nd);
    const double vx = sum_dx_d * inv_d;
    const double vy = sum_dy_d * inv_d;
    const double vz = sum_dz_d * inv_d;

    const double mx = sum_dx_sel * inv_sel;
    const double my = sum_dy_sel * inv_sel;
    const double mz = sum_dz_sel * inv_sel;

    // mean |dr - v|^2 = mean|dr|^2 - 2 (v dot mean(dr)) + |v|^2
    double v_dot_m = 0.0;
    double v2 = 0.0;
    if (use_x) { v_dot_m += vx*mx; v2 += vx*vx; }
    if (use_y) { v_dot_m += vy*my; v2 += vy*vy; }
    if (use_z) { v_dot_m += vz*mz; v2 += vz*vz; }
    const double msd = mean_dr2_sel - 2.0*v_dot_m + v2;
    return (msd > 0.0) ? msd : 0.0;
  }

private:
  int comp_mask_ = 7; // x|y|z
  bool initialized_ = false;
  std::int64_t origin_timestep_ = 0;
  std::vector<double> x0_, y0_, z0_;
};

} // namespace pilots
