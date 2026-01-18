#pragma once

#include <cstddef>
#include <functional>
#include <stdexcept>

#include "pilots/select/SelectionView.hpp"

#if defined(PILOTS_HAS_OPENMP) && PILOTS_HAS_OPENMP
  #include <omp.h>
#endif

namespace pilots::sdk {

// Parallel helpers for iterating selections.
//
// These are deliberately small and explicit: measures keep full control over
// algorithm structure; the SDK only provides safe OpenMP-friendly traversal
// primitives.

template <class F>
inline void for_each_atom(const pilots::SelectionView& sel, F&& fn) {
  const std::size_t n = sel.size();
  const auto* idx = sel.idx.data();
#if defined(PILOTS_HAS_OPENMP) && PILOTS_HAS_OPENMP
  #pragma omp parallel for schedule(static)
  for (std::int64_t ii = 0; ii < static_cast<std::int64_t>(n); ++ii) {
    const std::size_t a = idx[static_cast<std::size_t>(ii)];
    fn(a);
  }
#else
  for (std::size_t ii = 0; ii < n; ++ii) {
    fn(idx[ii]);
  }
#endif
}

// Sum reduction helper.
//
// If deterministic=true, we use a single-thread loop for strict reproducibility.
// Otherwise we use OpenMP reduction when available.

template <class F>
inline double reduce_sum_f64(const pilots::SelectionView& sel, F&& fn, bool deterministic = false) {
  const std::size_t n = sel.size();
  const auto* idx = sel.idx.data();
  if (deterministic) {
    double sum = 0.0;
    for (std::size_t ii = 0; ii < n; ++ii) sum += fn(idx[ii]);
    return sum;
  }

#if defined(PILOTS_HAS_OPENMP) && PILOTS_HAS_OPENMP
  double sum = 0.0;
  #pragma omp parallel for reduction(+:sum) schedule(static)
  for (std::int64_t ii = 0; ii < static_cast<std::int64_t>(n); ++ii) {
    sum += fn(idx[static_cast<std::size_t>(ii)]);
  }
  return sum;
#else
  double sum = 0.0;
  for (std::size_t ii = 0; ii < n; ++ii) sum += fn(idx[ii]);
  return sum;
#endif
}

} // namespace pilots::sdk
