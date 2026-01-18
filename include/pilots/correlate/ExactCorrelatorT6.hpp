#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#if PILOTS_HAS_OPENMP
#include <omp.h>
#endif

#include "pilots/correlate/CorrelationSeries.hpp"
#include "pilots/correlate/ICorrelatorT6.hpp"
#include "pilots/correlate/LagAxis.hpp"
#include "pilots/correlate/Tensor6Types.hpp"
#include "pilots/util/BinaryIO.hpp"

namespace pilots {

// Exact time-origin correlator for 6-channel tensor signals.
//
// PairOp must provide:
//   Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const;
// It may return a scalar observable by placing it in v[xx] and leaving other channels 0.

template <typename PairOp>
class ExactCorrelatorT6 final : public ICorrelatorT6 {
public:
  struct Options {
    double dt = 0.0;               // seconds per timestep
    LagAxis axis = LagAxis::Frame;
    double timebin_width = 0.0;    // seconds (only for axis=TimeBin)
    std::size_t window_frames = 0; // number of frames fed into this correlator
    std::size_t lag_stride = 1;    // only meaningful for axis=Frame
    std::size_t block_size = 200;  // 0 disables SEM

    // timebin post-processing
    std::uint64_t min_pairs_per_bin = 0; // 0 disables merging
    bool bin_merge = false;             // if true and min_pairs_per_bin>0, merge sparse bins
  };

  ExactCorrelatorT6(std::size_t nsel, Options opt, PairOp pair)
  : nsel_(nsel), opt_(opt), pair_(std::move(pair)) {
    if (!(opt_.dt > 0.0)) throw std::runtime_error("ExactCorrelatorT6: dt must be > 0");
    if (nsel_ == 0) throw std::runtime_error("ExactCorrelatorT6: nsel must be > 0");
    if (opt_.window_frames < 2) throw std::runtime_error("ExactCorrelatorT6: window_frames must be >= 2");
    if (opt_.lag_stride < 1) throw std::runtime_error("ExactCorrelatorT6: lag_stride must be >= 1");
    if (opt_.axis == LagAxis::TimeBin && !(opt_.timebin_width > 0.0)) {
      throw std::runtime_error("ExactCorrelatorT6: timebin_width must be > 0 for lag_axis=timebin");
    }
    if ((opt_.axis == LagAxis::Timestep || opt_.axis == LagAxis::TimeBin) && opt_.lag_stride != 1) {
      throw std::runtime_error("ExactCorrelatorT6: lag_stride is only supported for lag_axis=frame");
    }

    ring_.resize(opt_.window_frames);
    for (auto& s : ring_) {
      s.timestep = std::numeric_limits<std::int64_t>::min();
      s.xx.resize(nsel_);
      s.yy.resize(nsel_);
      s.zz.resize(nsel_);
      s.xy.resize(nsel_);
      s.xz.resize(nsel_);
      s.yz.resize(nsel_);
    }

    if (opt_.axis == LagAxis::Frame) {
      stats_vec_.assign(opt_.window_frames, LagStatsT6{});
    }
  }

  void start() override { started_ = true; }

  void push(std::int64_t timestep,
            std::span<const double> xx,
            std::span<const double> yy,
            std::span<const double> zz,
            std::span<const double> xy,
            std::span<const double> xz,
            std::span<const double> yz) override {
    if (!started_) throw std::runtime_error("ExactCorrelatorT6: start() must be called before push()");
    if (xx.size() != nsel_ || yy.size() != nsel_ || zz.size() != nsel_ || xy.size() != nsel_ || xz.size() != nsel_ || yz.size() != nsel_) {
      throw std::runtime_error("ExactCorrelatorT6: input vector sizes do not match nsel");
    }
    if (nsamples_ >= opt_.window_frames) {
      throw std::runtime_error("ExactCorrelatorT6: pushed more samples than window_frames");
    }

    const std::size_t cur_idx = nsamples_;
    T6Slot& cur = ring_[cur_idx];
    cur.timestep = timestep;

#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
    for (std::size_t p = 0; p < nsel_; ++p) {
      cur.xx[p] = xx[p];
      cur.yy[p] = yy[p];
      cur.zz[p] = zz[p];
      cur.xy[p] = xy[p];
      cur.xz[p] = xz[p];
      cur.yz[p] = yz[p];
    }

    // Enumerate all lags for this current frame
    for (std::size_t lag = opt_.lag_stride; lag <= cur_idx; lag += opt_.lag_stride) {
      const std::size_t org_idx = cur_idx - lag;
      const T6Slot& org = ring_[org_idx];
      if (org.timestep == std::numeric_limits<std::int64_t>::min()) continue;

      const std::int64_t dtstep_i = cur.timestep - org.timestep;
      if (dtstep_i <= 0) continue;

      const Tensor6 val = pair_(cur, org);

      if (opt_.axis == LagAxis::Frame) {
        LagStatsT6& st = stats_vec_[lag];
        update_stats_(st, val, static_cast<double>(dtstep_i));
      } else if (opt_.axis == LagAxis::Timestep) {
        LagStatsT6& st = stats_map_[dtstep_i];
        update_stats_(st, val, static_cast<double>(dtstep_i));
      } else { // TimeBin
        const double dt = static_cast<double>(dtstep_i) * opt_.dt;
        const std::int64_t bin = static_cast<std::int64_t>(std::floor(dt / opt_.timebin_width));
        LagStatsT6& st = stats_map_[bin];
        update_stats_(st, val, static_cast<double>(dtstep_i));
      }
    }

    nsamples_ += 1;
  }

  CorrelationSeriesT6 finalize() override {
    CorrelationSeriesT6 out;
    out.axis = opt_.axis;

    if (opt_.axis == LagAxis::Frame) {
      for (std::size_t lag = opt_.lag_stride; lag < opt_.window_frames; lag += opt_.lag_stride) {
        const LagStatsT6& st = stats_vec_[lag];
        if (st.count == 0) continue;
        const double mean_dtstep = st.sum_dtimestep / static_cast<double>(st.count);
        const double time = mean_dtstep * opt_.dt;
        append_entry_(out, static_cast<double>(lag), st, time, mean_dtstep);
      }
      return out;
    }

    // Map-based axes: timestep or timebin
    std::vector<std::pair<std::int64_t, LagStatsT6>> items;
    items.reserve(stats_map_.size());
    for (auto& kv : stats_map_) items.emplace_back(kv.first, kv.second);
    std::sort(items.begin(), items.end(), [](const auto& a, const auto& b){ return a.first < b.first; });

    if (opt_.axis == LagAxis::TimeBin && opt_.bin_merge && opt_.min_pairs_per_bin > 0) {
      items = merge_bins_(std::move(items));
    }

    for (const auto& kv : items) {
      const std::int64_t key_i = kv.first;
      const LagStatsT6& st = kv.second;
      if (st.count == 0) continue;

      if (opt_.axis == LagAxis::Timestep) {
        const double lag_key = static_cast<double>(key_i);
        const double mean_dtstep = st.sum_dtimestep / static_cast<double>(st.count);
        const double time = lag_key * opt_.dt;
        append_entry_(out, lag_key, st, time, mean_dtstep);
      } else {
        // timebin: key is bin index (informational); output uses mean time from data.
        const double mean_dtstep = st.sum_dtimestep / static_cast<double>(st.count);
        const double time = mean_dtstep * opt_.dt;
        append_entry_(out, time, st, time, mean_dtstep);
      }
    }

    return out;
  }

  CorrelationSeriesT6 snapshot() const override {
    return const_cast<ExactCorrelatorT6*>(this)->finalize();
  }

  void save_state(std::ostream& os) const override {
    util::write_magic(os, "PILOTSEXACTT6");
    util::BinaryWriter w(os);
    w.write_u32(1);
    w.write_u64(static_cast<std::uint64_t>(nsel_));

    w.write_f64(opt_.dt);
    w.write_i32(static_cast<std::int32_t>(opt_.axis));
    w.write_f64(opt_.timebin_width);
    w.write_u64(static_cast<std::uint64_t>(opt_.window_frames));
    w.write_u64(static_cast<std::uint64_t>(opt_.lag_stride));
    w.write_u64(static_cast<std::uint64_t>(opt_.block_size));
    w.write_u64(opt_.min_pairs_per_bin);
    w.write_u8(opt_.bin_merge ? 1 : 0);

    w.write_u8(started_ ? 1 : 0);
    w.write_u64(static_cast<std::uint64_t>(nsamples_));

    // samples
    for (std::size_t i = 0; i < nsamples_; ++i) {
      write_slot_(w, ring_[i]);
    }

    // stats (frame axis uses dense vector)
    w.write_u64(static_cast<std::uint64_t>(stats_vec_.size()));
    for (const auto& st : stats_vec_) write_lagstats_(w, st);

    // stats_map
    w.write_u64(static_cast<std::uint64_t>(stats_map_.size()));
    for (const auto& kv : stats_map_) {
      w.write_i64(kv.first);
      write_lagstats_(w, kv.second);
    }

    util::write_magic(os, "PILOTSEND");
  }

  void load_state(std::istream& is) override {
    util::require_magic(is, "PILOTSEXACTT6");
    util::BinaryReader r(is);
    const std::uint32_t ver = r.read_u32();
    if (ver != 1) throw std::runtime_error("ExactCorrelatorT6: unsupported checkpoint version");

    const std::uint64_t nsel_ck = r.read_u64();
    if (nsel_ck != static_cast<std::uint64_t>(nsel_)) throw std::runtime_error("ExactCorrelatorT6: nsel mismatch in checkpoint");

    auto near_eq = [](double a, double b) {
      const double d = (a > b) ? (a - b) : (b - a);
      return d <= 1e-12;
    };

    const double dt_ck = r.read_f64();
    const LagAxis axis_ck = static_cast<LagAxis>(r.read_i32());
    const double tb_ck = r.read_f64();
    const std::uint64_t win_ck = r.read_u64();
    const std::uint64_t stride_ck = r.read_u64();
    const std::uint64_t block_ck = r.read_u64();
    const std::uint64_t minpairs_ck = r.read_u64();
    const bool binmerge_ck = r.read_u8() != 0;

    if (!near_eq(dt_ck, opt_.dt)) throw std::runtime_error("ExactCorrelatorT6: dt mismatch in checkpoint");
    if (axis_ck != opt_.axis) throw std::runtime_error("ExactCorrelatorT6: axis mismatch in checkpoint");
    if (!near_eq(tb_ck, opt_.timebin_width)) throw std::runtime_error("ExactCorrelatorT6: timebin_width mismatch in checkpoint");
    if (win_ck != static_cast<std::uint64_t>(opt_.window_frames)) throw std::runtime_error("ExactCorrelatorT6: window_frames mismatch in checkpoint");
    if (stride_ck != static_cast<std::uint64_t>(opt_.lag_stride)) throw std::runtime_error("ExactCorrelatorT6: lag_stride mismatch in checkpoint");
    if (block_ck != static_cast<std::uint64_t>(opt_.block_size)) throw std::runtime_error("ExactCorrelatorT6: block_size mismatch in checkpoint");
    if (minpairs_ck != opt_.min_pairs_per_bin) throw std::runtime_error("ExactCorrelatorT6: min_pairs_per_bin mismatch in checkpoint");
    if (binmerge_ck != opt_.bin_merge) throw std::runtime_error("ExactCorrelatorT6: bin_merge mismatch in checkpoint");

    started_ = (r.read_u8() != 0);
    const std::uint64_t ns_ck = r.read_u64();
    if (ns_ck > static_cast<std::uint64_t>(opt_.window_frames)) throw std::runtime_error("ExactCorrelatorT6: nsamples exceeds window_frames in checkpoint");
    nsamples_ = static_cast<std::size_t>(ns_ck);

    // reset ring
    for (auto& s : ring_) {
      s.timestep = std::numeric_limits<std::int64_t>::min();
    }
    for (std::size_t i = 0; i < nsamples_; ++i) {
      read_slot_(r, ring_[i]);
    }

    // stats_vec
    const std::uint64_t sv_sz = r.read_u64();
    if (sv_sz != static_cast<std::uint64_t>(stats_vec_.size())) {
      throw std::runtime_error("ExactCorrelatorT6: stats_vec size mismatch in checkpoint");
    }
    for (std::size_t i = 0; i < stats_vec_.size(); ++i) {
      read_lagstats_(r, stats_vec_[i]);
    }

    // stats_map
    const std::uint64_t sm_sz = r.read_u64();
    stats_map_.clear();
    stats_map_.reserve(static_cast<std::size_t>(sm_sz) * 2);
    for (std::uint64_t i = 0; i < sm_sz; ++i) {
      const std::int64_t key = r.read_i64();
      LagStatsT6 st;
      read_lagstats_(r, st);
      stats_map_.emplace(key, st);
    }

    util::require_magic(is, "PILOTSEND");
  }

private:
  std::size_t nsel_ = 0;
  Options opt_;
  PairOp pair_;
  bool started_ = false;
  std::size_t nsamples_ = 0;

  std::vector<T6Slot> ring_;
  std::vector<LagStatsT6> stats_vec_;
  std::unordered_map<std::int64_t, LagStatsT6> stats_map_;

  void update_stats_(LagStatsT6& st, const Tensor6& val, double dtstep) {
    for (int c = 0; c < 6; ++c) {
      st.sum_val[static_cast<std::size_t>(c)] += val.v[static_cast<std::size_t>(c)];
    }
    st.sum_dtimestep += dtstep;
    st.count += 1;

    if (opt_.block_size == 0) return;

    for (int c = 0; c < 6; ++c) {
      st.block_cur_sum[static_cast<std::size_t>(c)] += val.v[static_cast<std::size_t>(c)];
    }
    st.block_cur_count += 1;
    if (st.block_cur_count >= opt_.block_size) {
      const double denom = static_cast<double>(st.block_cur_count);
      for (int c = 0; c < 6; ++c) {
        const double bm = st.block_cur_sum[static_cast<std::size_t>(c)] / denom;
        st.block_mean_sum[static_cast<std::size_t>(c)] += bm;
        st.block_mean_sumsq[static_cast<std::size_t>(c)] += bm * bm;
        st.block_cur_sum[static_cast<std::size_t>(c)] = 0.0;
      }
      st.block_count += 1;
      st.block_cur_count = 0;
    }
  }

  static void merge_stats_into_(LagStatsT6& dst, const LagStatsT6& src) {
    for (int c = 0; c < 6; ++c) {
      dst.sum_val[static_cast<std::size_t>(c)] += src.sum_val[static_cast<std::size_t>(c)];
      dst.block_mean_sum[static_cast<std::size_t>(c)] += src.block_mean_sum[static_cast<std::size_t>(c)];
      dst.block_mean_sumsq[static_cast<std::size_t>(c)] += src.block_mean_sumsq[static_cast<std::size_t>(c)];
      // NOTE: we intentionally do not merge partial block_cur_* to avoid cross-bin mixing.
    }
    dst.sum_dtimestep += src.sum_dtimestep;
    dst.count += src.count;
    dst.block_count += src.block_count;
  }

  std::vector<std::pair<std::int64_t, LagStatsT6>> merge_bins_(std::vector<std::pair<std::int64_t, LagStatsT6>> items) const {
    std::vector<std::pair<std::int64_t, LagStatsT6>> out;
    out.reserve(items.size());
    const std::uint64_t thr = opt_.min_pairs_per_bin;
    std::size_t i = 0;
    while (i < items.size()) {
      std::int64_t key = items[i].first;
      LagStatsT6 st = items[i].second;
      while (st.count < thr && (i + 1) < items.size()) {
        ++i;
        // Keep the first key as an identifier; output time uses mean_dtstep so key is informational.
        merge_stats_into_(st, items[i].second);
      }
      out.emplace_back(key, st);
      ++i;
    }

    // If the last bin is still under threshold, merge it backward.
    if (out.size() >= 2 && out.back().second.count < thr) {
      auto last = out.back();
      out.pop_back();
      merge_stats_into_(out.back().second, last.second);
    }
    return out;
  }

  void append_entry_(CorrelationSeriesT6& out, double lag_key, const LagStatsT6& st, double time, double mean_dtstep) const {
    const double invc = 1.0 / static_cast<double>(st.count);
    const double mean_xx = st.sum_val[T6_XX] * invc;
    const double mean_yy = st.sum_val[T6_YY] * invc;
    const double mean_zz = st.sum_val[T6_ZZ] * invc;
    const double mean_xy = st.sum_val[T6_XY] * invc;
    const double mean_xz = st.sum_val[T6_XZ] * invc;
    const double mean_yz = st.sum_val[T6_YZ] * invc;

    std::array<double,6> sem{0.0,0.0,0.0,0.0,0.0,0.0};
    std::uint64_t nb = 0;
    if (opt_.block_size > 0) {
      nb = st.block_count;
      if (nb >= 2) {
        for (int c = 0; c < 6; ++c) {
          const double mb = st.block_mean_sum[static_cast<std::size_t>(c)] / static_cast<double>(nb);
          const double m2b = st.block_mean_sumsq[static_cast<std::size_t>(c)] / static_cast<double>(nb);
          double var = m2b - mb * mb;
          if (var < 0.0) var = 0.0;
          var *= static_cast<double>(nb) / static_cast<double>(nb - 1);
          sem[static_cast<std::size_t>(c)] = std::sqrt(var / static_cast<double>(nb));
        }
      }
    }

    out.lag.push_back(lag_key);
    out.time.push_back(time);
    out.value_xx.push_back(mean_xx);
    out.value_yy.push_back(mean_yy);
    out.value_zz.push_back(mean_zz);
    out.value_xy.push_back(mean_xy);
    out.value_xz.push_back(mean_xz);
    out.value_yz.push_back(mean_yz);
    out.sem_xx.push_back(sem[T6_XX]);
    out.sem_yy.push_back(sem[T6_YY]);
    out.sem_zz.push_back(sem[T6_ZZ]);
    out.sem_xy.push_back(sem[T6_XY]);
    out.sem_xz.push_back(sem[T6_XZ]);
    out.sem_yz.push_back(sem[T6_YZ]);
    out.count_pairs.push_back(st.count);
    out.n_blocks.push_back(nb);
    out.mean_dtimestep.push_back(mean_dtstep);
  }

  void write_slot_(util::BinaryWriter& w, const T6Slot& s) const {
    if (s.xx.size() != nsel_) throw std::runtime_error("ExactCorrelatorT6: slot size mismatch");
    w.write_i64(s.timestep);
    w.write_bytes(s.xx.data(), sizeof(double) * s.xx.size());
    w.write_bytes(s.yy.data(), sizeof(double) * s.yy.size());
    w.write_bytes(s.zz.data(), sizeof(double) * s.zz.size());
    w.write_bytes(s.xy.data(), sizeof(double) * s.xy.size());
    w.write_bytes(s.xz.data(), sizeof(double) * s.xz.size());
    w.write_bytes(s.yz.data(), sizeof(double) * s.yz.size());
  }

  void read_slot_(util::BinaryReader& r, T6Slot& s) const {
    s.timestep = r.read_i64();
    if (s.xx.size() != nsel_) throw std::runtime_error("ExactCorrelatorT6: slot size mismatch on load");
    r.read_bytes(s.xx.data(), sizeof(double) * s.xx.size());
    r.read_bytes(s.yy.data(), sizeof(double) * s.yy.size());
    r.read_bytes(s.zz.data(), sizeof(double) * s.zz.size());
    r.read_bytes(s.xy.data(), sizeof(double) * s.xy.size());
    r.read_bytes(s.xz.data(), sizeof(double) * s.xz.size());
    r.read_bytes(s.yz.data(), sizeof(double) * s.yz.size());
  }

  static void write_lagstats_(util::BinaryWriter& w, const LagStatsT6& st) {
    for (int c = 0; c < 6; ++c) w.write_f64(st.sum_val[static_cast<std::size_t>(c)]);
    w.write_f64(st.sum_dtimestep);
    w.write_u64(st.count);
    for (int c = 0; c < 6; ++c) w.write_f64(st.block_cur_sum[static_cast<std::size_t>(c)]);
    w.write_u64(st.block_cur_count);
    for (int c = 0; c < 6; ++c) w.write_f64(st.block_mean_sum[static_cast<std::size_t>(c)]);
    for (int c = 0; c < 6; ++c) w.write_f64(st.block_mean_sumsq[static_cast<std::size_t>(c)]);
    w.write_u64(st.block_count);
  }

  static void read_lagstats_(util::BinaryReader& r, LagStatsT6& st) {
    for (int c = 0; c < 6; ++c) st.sum_val[static_cast<std::size_t>(c)] = r.read_f64();
    st.sum_dtimestep = r.read_f64();
    st.count = r.read_u64();
    for (int c = 0; c < 6; ++c) st.block_cur_sum[static_cast<std::size_t>(c)] = r.read_f64();
    st.block_cur_count = r.read_u64();
    for (int c = 0; c < 6; ++c) st.block_mean_sum[static_cast<std::size_t>(c)] = r.read_f64();
    for (int c = 0; c < 6; ++c) st.block_mean_sumsq[static_cast<std::size_t>(c)] = r.read_f64();
    st.block_count = r.read_u64();
  }
};

} // namespace pilots
