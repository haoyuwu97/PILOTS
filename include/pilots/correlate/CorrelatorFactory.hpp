#pragma once

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

#include "pilots/config/IniConfig.hpp"
#include "pilots/correlate/CorrelatorSpec.hpp"
#include "pilots/correlate/ExactCorrelatorT6.hpp"
#include "pilots/correlate/MultiTauCorrelatorT6.hpp"
#include "pilots/correlate/FFTCorrelatorStubT6.hpp"
#include "pilots/correlate/T6AvgOp.hpp"

namespace pilots {

inline CorrelatorSpec parse_correlator_spec(const IniConfig& cfg, const std::string& section, double dt) {
  CorrelatorSpec s;
  s.dt = dt;

  std::string t = cfg.get_string(section, "correlator", std::optional<std::string>("exact"));
  for (auto& c : t) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  s.type = t;

  if (cfg.has_key(section, "lag_axis")) {
    s.axis = parse_lag_axis(cfg.get_string(section, "lag_axis"));
  }

  if (s.axis == LagAxis::TimeBin) {
    s.timebin_width = cfg.get_double(section, "timebin_width", std::optional<double>(dt));
    if (!(s.timebin_width > 0.0)) {
      throw std::runtime_error("timebin_width must be > 0 for lag_axis=timebin");
    }

    s.min_pairs_per_bin = static_cast<std::uint64_t>(cfg.get_int64(section, "min_pairs_per_bin", std::optional<std::int64_t>(0)));
    const bool has_bin_merge = cfg.has_key(section, "bin_merge");
    if (has_bin_merge) {
      s.bin_merge = cfg.get_bool(section, "bin_merge", std::optional<bool>(false));
    } else {
      s.bin_merge = (s.min_pairs_per_bin > 0);
    }
  }

  s.lag_stride = cfg.get_size(section, "lag_stride", std::optional<std::size_t>(1));
  s.block_size = cfg.get_size(section, "block_size", std::optional<std::size_t>(200));
  s.mt_channels = cfg.get_size(section, "mt_channels", std::optional<std::size_t>(16));
  s.mt_levels = cfg.get_size(section, "mt_levels", std::optional<std::size_t>(10));

  return s;
}

template <typename PairOp, typename AvgOp>
std::unique_ptr<ICorrelatorT6> make_correlator_t6(std::size_t nsel,
                                                  std::size_t window_frames,
                                                  const CorrelatorSpec& spec,
                                                  PairOp pair,
                                                  AvgOp avg) {
  if (spec.type == "exact") {
    typename ExactCorrelatorT6<PairOp>::Options opt;
    opt.dt = spec.dt;
    opt.axis = spec.axis;
    opt.timebin_width = spec.timebin_width;
    opt.window_frames = window_frames;
    opt.lag_stride = spec.lag_stride;
    opt.block_size = spec.block_size;
    opt.min_pairs_per_bin = spec.min_pairs_per_bin;
    opt.bin_merge = spec.bin_merge;
    return std::make_unique<ExactCorrelatorT6<PairOp>>(nsel, opt, std::move(pair));
  }
  if (spec.type == "multitau") {
    if (spec.axis != LagAxis::Frame) {
      throw std::runtime_error("correlator=multitau currently supports lag_axis=frame only. Use correlator=exact with lag_axis=timestep/timebin for irregular sampling.");
    }
    typename MultiTauCorrelatorT6<PairOp, AvgOp>::Options opt;
    opt.dt = spec.dt;
    opt.channels = spec.mt_channels;
    opt.levels = spec.mt_levels;
    opt.block_size = spec.block_size;
    return std::make_unique<MultiTauCorrelatorT6<PairOp, AvgOp>>(nsel, opt, std::move(pair), std::move(avg));
  }
  if (spec.type == "fft") {
    typename FFTCorrelatorStubT6<PairOp>::Options opt;
    opt.dt = spec.dt;
    return std::make_unique<FFTCorrelatorStubT6<PairOp>>(nsel, opt, std::move(pair));
  }
  throw std::runtime_error("unsupported correlator type: " + spec.type);
}

// ---- Runtime (non-template) factory ----
//
// Motivation:
//   Future measures should not need to mention template parameters for correlators.
//   They can provide PairOp/AvgOp behavior through std::function (type erasure), and
//   the factory instantiates correlators with a fixed runtime op type.
//
// Pair function: two samples -> Tensor6 observable
using PairFnT6 = std::function<Tensor6(const T6Slot& cur, const T6Slot& org)>;
// Avg function (for multi-tau): two samples -> coarse sample
using AvgFnT6  = std::function<void(const T6Slot& a, const T6Slot& b, T6Slot& out)>;

struct PairOpRuntimeT6 {
  PairFnT6 fn;
  Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const {
    if (!fn) throw std::runtime_error("PairOpRuntimeT6: pair function not set");
    return fn(cur, org);
  }
};

struct AvgOpRuntimeT6 {
  AvgFnT6 fn;
  void operator()(const T6Slot& a, const T6Slot& b, T6Slot& out) const {
    if (!fn) throw std::runtime_error("AvgOpRuntimeT6: avg function not set");
    fn(a, b, out);
  }
};

// Non-template correlator factory. Measures provide `pair_fn` (always) and optionally
// `avg_fn` (required only when correlator=multitau). The returned correlator is fully
// configured by `spec`.
inline std::unique_ptr<ICorrelatorT6> make_correlator_t6_runtime(std::size_t nsel,
                                                                std::size_t window_frames,
                                                                const CorrelatorSpec& spec,
                                                                PairFnT6 pair_fn,
                                                                AvgFnT6 avg_fn = {}) {
  if (!pair_fn) {
    throw std::runtime_error("make_correlator_t6_runtime: pair_fn must be provided");
  }
  if (spec.type == "multitau" && !avg_fn) {
    throw std::runtime_error("make_correlator_t6_runtime: avg_fn must be provided for correlator=multitau");
  }

  PairOpRuntimeT6 pair{std::move(pair_fn)};
  AvgOpRuntimeT6 avg{std::move(avg_fn)};
  return make_correlator_t6<PairOpRuntimeT6, AvgOpRuntimeT6>(nsel, window_frames, spec, std::move(pair), std::move(avg));
}

// Convenience overload: measures provide a PairOp object (a functor with
//   Tensor6 operator()(const T6Slot&, const T6Slot&) const
// and the factory performs type-erasure automatically.
//
// For correlator=multitau, a default AvgOp (element-wise T6AvgOp) is used unless
// the caller overrides it.
template <typename PairOp>
inline std::unique_ptr<ICorrelatorT6> make_correlator_t6_runtime_auto(std::size_t nsel,
                                                                      std::size_t window_frames,
                                                                      const CorrelatorSpec& spec,
                                                                      PairOp pair_op,
                                                                      AvgFnT6 avg_override = {}) {
  auto pair_ptr = std::make_shared<PairOp>(std::move(pair_op));
  PairFnT6 pair_fn = [pair_ptr](const T6Slot& cur, const T6Slot& org) -> Tensor6 {
    return (*pair_ptr)(cur, org);
  };

  AvgFnT6 avg_fn;
  if (spec.type == "multitau") {
    if (avg_override) {
      avg_fn = std::move(avg_override);
    } else {
      auto avg_ptr = std::make_shared<T6AvgOp>();
      avg_fn = [avg_ptr](const T6Slot& a, const T6Slot& b, T6Slot& out) {
        (*avg_ptr)(a, b, out);
      };
    }
  }

  return make_correlator_t6_runtime(nsel, window_frames, spec, std::move(pair_fn), std::move(avg_fn));
}

} // namespace pilots
