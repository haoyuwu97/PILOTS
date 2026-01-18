#pragma once

#include <cstddef>
#include <stdexcept>

#if PILOTS_HAS_OPENMP
#include <omp.h>
#endif

#include "pilots/correlate/Tensor6Types.hpp"

namespace pilots {

// Element-wise average of two T6Slot samples into an output slot.
// `out` must already have all channels resized to the same length.
struct T6AvgOp {
  void operator()(const T6Slot& a, const T6Slot& b, T6Slot& out) const {
    const std::size_t n = a.xx.size();
    if (b.xx.size() != n || out.xx.size() != n ||
        a.yy.size() != n || b.yy.size() != n || out.yy.size() != n ||
        a.zz.size() != n || b.zz.size() != n || out.zz.size() != n ||
        a.xy.size() != n || b.xy.size() != n || out.xy.size() != n ||
        a.xz.size() != n || b.xz.size() != n || out.xz.size() != n ||
        a.yz.size() != n || b.yz.size() != n || out.yz.size() != n) {
      throw std::runtime_error("T6AvgOp: slot sizes mismatch");
    }

    // integer average to avoid overflow
    out.timestep = (a.timestep / 2) + (b.timestep / 2);

#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
    for (std::size_t p = 0; p < n; ++p) {
      out.xx[p] = 0.5 * (a.xx[p] + b.xx[p]);
      out.yy[p] = 0.5 * (a.yy[p] + b.yy[p]);
      out.zz[p] = 0.5 * (a.zz[p] + b.zz[p]);
      out.xy[p] = 0.5 * (a.xy[p] + b.xy[p]);
      out.xz[p] = 0.5 * (a.xz[p] + b.xz[p]);
      out.yz[p] = 0.5 * (a.yz[p] + b.yz[p]);
    }
  }
};

} // namespace pilots
