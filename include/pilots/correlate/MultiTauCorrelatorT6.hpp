#pragma once

#include <algorithm>
#include <array>
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

// Multi-tau time-origin correlator for T6Slot samples.
//
// PairOp must provide:
//   Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const;
// AvgOp must provide:
//   void operator()(const T6Slot& a, const T6Slot& b, T6Slot& out) const;
// where `out` already has all channels resized to nsel.

template <typename PairOp, typename AvgOp>
class MultiTauCorrelatorT6 final : public ICorrelatorT6 {
public:
  struct Options {
    double dt = 0.0;          // seconds per timestep
    std::size_t channels = 16; // must be even
    std::size_t levels = 10;   // >=1
    std::size_t block_size = 200; // 0 disables SEM
  };

  MultiTauCorrelatorT6(std::size_t nsel, Options opt, PairOp pair, AvgOp avg)
  : nsel_(nsel), opt_(opt), pair_(std::move(pair)), avg_(std::move(avg)) {
    if (!(opt_.dt > 0.0)) throw std::runtime_error("MultiTauCorrelatorT6: dt must be > 0");
    if (nsel_ == 0) throw std::runtime_error("MultiTauCorrelatorT6: nsel must be > 0");
    if (opt_.channels < 2 || (opt_.channels % 2) != 0) throw std::runtime_error("MultiTauCorrelatorT6: channels must be an even integer >= 2");
    if (opt_.levels < 1) throw std::runtime_error("MultiTauCorrelatorT6: levels must be >= 1");
    build_lag_schedule_();
    allocate_stats_();
  }

  void start() override {
    started_ = true;
    levels_.resize(opt_.levels);
    for (std::size_t lev = 0; lev < opt_.levels; ++lev) {
      Level& L = levels_[lev];
      L.m = opt_.channels;
      L.buf.resize(L.m + 1);
      for (auto& s : L.buf) {
        s.timestep = std::numeric_limits<std::int64_t>::min();
        resize_slot_(s);
      }
      L.tmp_avg.timestep = std::numeric_limits<std::int64_t>::min();
      resize_slot_(L.tmp_avg);

      L.write_pos = 0;
      L.filled = 0;
      L.pending_idx = -1;
    }
  }

  void push(std::int64_t timestep,
            std::span<const double> xx,
            std::span<const double> yy,
            std::span<const double> zz,
            std::span<const double> xy,
            std::span<const double> xz,
            std::span<const double> yz) override {
    if (!started_) throw std::runtime_error("MultiTauCorrelatorT6: start() must be called before push()");
    if (xx.size() != nsel_ || yy.size() != nsel_ || zz.size() != nsel_ || xy.size() != nsel_ || xz.size() != nsel_ || yz.size() != nsel_) {
      throw std::runtime_error("MultiTauCorrelatorT6: input vector sizes do not match nsel");
    }

    Level& L0 = levels_[0];
    const std::size_t cur_idx = L0.write_pos;
    T6Slot& cur = L0.buf[cur_idx];
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

    insert_and_accumulate_(0, static_cast<int>(cur_idx));
  }

  CorrelationSeriesT6 finalize() override {
    CorrelationSeriesT6 out;
    out.axis = LagAxis::Frame;

    for (std::size_t k = 0; k < lags_.size(); ++k) {
      const std::uint64_t c = count_[k];
      if (c == 0) continue;
      const double invc = 1.0 / static_cast<double>(c);
      const double mean_dtstep = sum_dtimestep_[k] * invc;
      const double t = mean_dtstep * opt_.dt;

      std::array<double,6> mean{0,0,0,0,0,0};
      for (int ch = 0; ch < 6; ++ch) {
        mean[static_cast<std::size_t>(ch)] = sum_val_[k][static_cast<std::size_t>(ch)] * invc;
      }

      std::array<double,6> sem{0,0,0,0,0,0};
      std::uint64_t nb = 0;
      if (opt_.block_size > 0) {
        nb = block_count_[k];
        if (nb >= 2) {
          const double invnb = 1.0 / static_cast<double>(nb);
          for (int ch = 0; ch < 6; ++ch) {
            const double mb = block_mean_sum_[k][static_cast<std::size_t>(ch)] * invnb;
            const double m2b = block_mean_sumsq_[k][static_cast<std::size_t>(ch)] * invnb;
            double var = m2b - mb * mb;
            if (var < 0.0) var = 0.0;
            var *= static_cast<double>(nb) / static_cast<double>(nb - 1);
            sem[static_cast<std::size_t>(ch)] = std::sqrt(var / static_cast<double>(nb));
          }
        }
      }

      out.lag.push_back(static_cast<double>(lags_[k]));
      out.time.push_back(t);
      out.value_xx.push_back(mean[T6_XX]);
      out.value_yy.push_back(mean[T6_YY]);
      out.value_zz.push_back(mean[T6_ZZ]);
      out.value_xy.push_back(mean[T6_XY]);
      out.value_xz.push_back(mean[T6_XZ]);
      out.value_yz.push_back(mean[T6_YZ]);
      out.sem_xx.push_back(sem[T6_XX]);
      out.sem_yy.push_back(sem[T6_YY]);
      out.sem_zz.push_back(sem[T6_ZZ]);
      out.sem_xy.push_back(sem[T6_XY]);
      out.sem_xz.push_back(sem[T6_XZ]);
      out.sem_yz.push_back(sem[T6_YZ]);
      out.count_pairs.push_back(c);
      out.n_blocks.push_back(nb);
      out.mean_dtimestep.push_back(mean_dtstep);
    }

    return out;
  }

  CorrelationSeriesT6 snapshot() const override {
    return const_cast<MultiTauCorrelatorT6*>(this)->finalize();
  }

  void save_state(std::ostream& os) const override {
    util::write_magic(os, "PILOTSMULTITAUT6");
    util::BinaryWriter w(os);
    w.write_u32(1);
    w.write_u64(static_cast<std::uint64_t>(nsel_));

    w.write_f64(opt_.dt);
    w.write_u64(static_cast<std::uint64_t>(opt_.channels));
    w.write_u64(static_cast<std::uint64_t>(opt_.levels));
    w.write_u64(static_cast<std::uint64_t>(opt_.block_size));

    w.write_u8(started_ ? 1 : 0);

    // lag schedule (sanity)
    w.write_u64(static_cast<std::uint64_t>(lags_.size()));
    for (std::size_t v : lags_) w.write_u64(static_cast<std::uint64_t>(v));

    // stats arrays
    const std::size_t K = lags_.size();
    for (std::size_t k = 0; k < K; ++k) {
      for (int ch = 0; ch < 6; ++ch) w.write_f64(sum_val_[k][static_cast<std::size_t>(ch)]);
      w.write_f64(sum_dtimestep_[k]);
      w.write_u64(count_[k]);

      for (int ch = 0; ch < 6; ++ch) w.write_f64(block_cur_sum_[k][static_cast<std::size_t>(ch)]);
      w.write_u64(block_cur_count_[k]);
      for (int ch = 0; ch < 6; ++ch) w.write_f64(block_mean_sum_[k][static_cast<std::size_t>(ch)]);
      for (int ch = 0; ch < 6; ++ch) w.write_f64(block_mean_sumsq_[k][static_cast<std::size_t>(ch)]);
      w.write_u64(block_count_[k]);
    }

    // levels (buffers)
    w.write_u64(static_cast<std::uint64_t>(levels_.size()));
    for (const auto& L : levels_) {
      w.write_u64(static_cast<std::uint64_t>(L.m));
      w.write_u64(static_cast<std::uint64_t>(L.write_pos));
      w.write_u64(static_cast<std::uint64_t>(L.filled));
      w.write_i32(static_cast<std::int32_t>(L.pending_idx));

      w.write_u64(static_cast<std::uint64_t>(L.buf.size()));
      for (const auto& s : L.buf) {
        write_slot_(w, s);
      }
      write_slot_(w, L.tmp_avg);
    }

    util::write_magic(os, "PILOTSEND");
  }

  void load_state(std::istream& is) override {
    util::require_magic(is, "PILOTSMULTITAUT6");
    util::BinaryReader r(is);
    const std::uint32_t ver = r.read_u32();
    if (ver != 1) throw std::runtime_error("MultiTauCorrelatorT6: unsupported checkpoint version");

    const std::uint64_t nsel_ck = r.read_u64();
    if (nsel_ck != static_cast<std::uint64_t>(nsel_)) throw std::runtime_error("MultiTauCorrelatorT6: nsel mismatch in checkpoint");

    auto near_eq = [](double a, double b) {
      const double d = (a > b) ? (a - b) : (b - a);
      return d <= 1e-12;
    };

    const double dt_ck = r.read_f64();
    const std::uint64_t ch_ck = r.read_u64();
    const std::uint64_t lv_ck = r.read_u64();
    const std::uint64_t bs_ck = r.read_u64();

    if (!near_eq(dt_ck, opt_.dt)) throw std::runtime_error("MultiTauCorrelatorT6: dt mismatch in checkpoint");
    if (ch_ck != static_cast<std::uint64_t>(opt_.channels)) throw std::runtime_error("MultiTauCorrelatorT6: channels mismatch in checkpoint");
    if (lv_ck != static_cast<std::uint64_t>(opt_.levels)) throw std::runtime_error("MultiTauCorrelatorT6: levels mismatch in checkpoint");
    if (bs_ck != static_cast<std::uint64_t>(opt_.block_size)) throw std::runtime_error("MultiTauCorrelatorT6: block_size mismatch in checkpoint");

    started_ = (r.read_u8() != 0);

    // lag schedule
    const std::uint64_t lsz = r.read_u64();
    if (lsz != static_cast<std::uint64_t>(lags_.size())) throw std::runtime_error("MultiTauCorrelatorT6: lag schedule size mismatch in checkpoint");
    for (std::size_t i = 0; i < lags_.size(); ++i) {
      const std::uint64_t v = r.read_u64();
      if (v != static_cast<std::uint64_t>(lags_[i])) throw std::runtime_error("MultiTauCorrelatorT6: lag schedule mismatch in checkpoint");
    }

    // stats arrays
    const std::size_t K = lags_.size();
    for (std::size_t k = 0; k < K; ++k) {
      for (int ch = 0; ch < 6; ++ch) sum_val_[k][static_cast<std::size_t>(ch)] = r.read_f64();
      sum_dtimestep_[k] = r.read_f64();
      count_[k] = r.read_u64();

      for (int ch = 0; ch < 6; ++ch) block_cur_sum_[k][static_cast<std::size_t>(ch)] = r.read_f64();
      block_cur_count_[k] = r.read_u64();
      for (int ch = 0; ch < 6; ++ch) block_mean_sum_[k][static_cast<std::size_t>(ch)] = r.read_f64();
      for (int ch = 0; ch < 6; ++ch) block_mean_sumsq_[k][static_cast<std::size_t>(ch)] = r.read_f64();
      block_count_[k] = r.read_u64();
    }

    // levels
    const std::uint64_t nlev = r.read_u64();
    if (nlev != static_cast<std::uint64_t>(levels_.size())) throw std::runtime_error("MultiTauCorrelatorT6: levels_.size mismatch in checkpoint (start() must have been called)");

    for (std::size_t lev = 0; lev < levels_.size(); ++lev) {
      Level& L = levels_[lev];
      const std::uint64_t m_ck = r.read_u64();
      const std::uint64_t wp_ck = r.read_u64();
      const std::uint64_t filled_ck = r.read_u64();
      const int pending_ck = r.read_i32();

      if (m_ck != static_cast<std::uint64_t>(L.m)) throw std::runtime_error("MultiTauCorrelatorT6: level.m mismatch in checkpoint");
      if (wp_ck >= static_cast<std::uint64_t>(L.buf.size())) throw std::runtime_error("MultiTauCorrelatorT6: level.write_pos out of range in checkpoint");
      if (filled_ck > static_cast<std::uint64_t>(L.buf.size())) throw std::runtime_error("MultiTauCorrelatorT6: level.filled out of range in checkpoint");

      L.write_pos = static_cast<std::size_t>(wp_ck);
      L.filled = static_cast<std::size_t>(filled_ck);
      L.pending_idx = pending_ck;

      const std::uint64_t bsz = r.read_u64();
      if (bsz != static_cast<std::uint64_t>(L.buf.size())) throw std::runtime_error("MultiTauCorrelatorT6: level.buf size mismatch in checkpoint");
      for (auto& s : L.buf) {
        read_slot_(r, s);
      }
      read_slot_(r, L.tmp_avg);
    }

    util::require_magic(is, "PILOTSEND");
  }

private:
  struct Level {
    std::size_t m = 0;
    std::vector<T6Slot> buf;
    T6Slot tmp_avg;
    std::size_t write_pos = 0;
    std::size_t filled = 0;
    int pending_idx = -1;
  };

  std::size_t nsel_ = 0;
  Options opt_;
  PairOp pair_;
  AvgOp avg_;
  bool started_ = false;

  std::vector<Level> levels_;

  std::vector<std::size_t> lags_;
  std::unordered_map<std::size_t, std::size_t> lag2idx_;

  std::vector<std::array<double,6>> sum_val_;
  std::vector<double> sum_dtimestep_;
  std::vector<std::uint64_t> count_;

  std::vector<std::array<double,6>> block_cur_sum_;
  std::vector<std::uint64_t> block_cur_count_;
  std::vector<std::array<double,6>> block_mean_sum_;
  std::vector<std::array<double,6>> block_mean_sumsq_;
  std::vector<std::uint64_t> block_count_;

  void resize_slot_(T6Slot& s) {
    s.xx.resize(nsel_);
    s.yy.resize(nsel_);
    s.zz.resize(nsel_);
    s.xy.resize(nsel_);
    s.xz.resize(nsel_);
    s.yz.resize(nsel_);
  }

  void build_lag_schedule_() {
    const std::size_t m = opt_.channels;
    const std::size_t half = m / 2;

    std::vector<std::size_t> lags;
    lags.reserve(opt_.levels * m);

    for (std::size_t lev = 0; lev < opt_.levels; ++lev) {
      const std::size_t scale = (std::size_t(1) << lev);
      if (lev == 0) {
        for (std::size_t d = 1; d <= m; ++d) lags.push_back(d * scale);
      } else {
        for (std::size_t d = half + 1; d <= m; ++d) lags.push_back(d * scale);
      }
    }

    std::sort(lags.begin(), lags.end());
    lags.erase(std::unique(lags.begin(), lags.end()), lags.end());
    lags_ = std::move(lags);

    lag2idx_.reserve(lags_.size() * 2);
    for (std::size_t i = 0; i < lags_.size(); ++i) lag2idx_[lags_[i]] = i;
  }

  void allocate_stats_() {
    const std::size_t K = lags_.size();
    sum_val_.assign(K, std::array<double,6>{0,0,0,0,0,0});
    sum_dtimestep_.assign(K, 0.0);
    count_.assign(K, 0);

    block_cur_sum_.assign(K, std::array<double,6>{0,0,0,0,0,0});
    block_cur_count_.assign(K, 0);
    block_mean_sum_.assign(K, std::array<double,6>{0,0,0,0,0,0});
    block_mean_sumsq_.assign(K, std::array<double,6>{0,0,0,0,0,0});
    block_count_.assign(K, 0);
  }

  void update_stats_(std::size_t k, const Tensor6& val, double dtstep) {
    for (int ch = 0; ch < 6; ++ch) {
      sum_val_[k][static_cast<std::size_t>(ch)] += val.v[static_cast<std::size_t>(ch)];
    }
    sum_dtimestep_[k] += dtstep;
    count_[k] += 1;

    if (opt_.block_size == 0) return;

    for (int ch = 0; ch < 6; ++ch) {
      block_cur_sum_[k][static_cast<std::size_t>(ch)] += val.v[static_cast<std::size_t>(ch)];
    }
    block_cur_count_[k] += 1;
    if (block_cur_count_[k] >= opt_.block_size) {
      const double denom = static_cast<double>(block_cur_count_[k]);
      for (int ch = 0; ch < 6; ++ch) {
        const double bm = block_cur_sum_[k][static_cast<std::size_t>(ch)] / denom;
        block_mean_sum_[k][static_cast<std::size_t>(ch)] += bm;
        block_mean_sumsq_[k][static_cast<std::size_t>(ch)] += bm * bm;
        block_cur_sum_[k][static_cast<std::size_t>(ch)] = 0.0;
      }
      block_count_[k] += 1;
      block_cur_count_[k] = 0;
    }
  }

  void insert_and_accumulate_(std::size_t lev, int cur_idx) {
    Level& L = levels_[lev];
    const std::size_t m = L.m;
    const std::size_t half = m / 2;

    L.write_pos = (static_cast<std::size_t>(cur_idx) + 1) % (m + 1);
    if (L.filled < (m + 1)) L.filled += 1;

    const T6Slot& cur = L.buf[static_cast<std::size_t>(cur_idx)];

    std::size_t d_min = 1;
    std::size_t d_max = m;
    if (lev > 0) d_min = half + 1;

    for (std::size_t d = d_min; d <= d_max; ++d) {
      if (L.filled <= d) continue;
      const std::size_t org_idx = (static_cast<std::size_t>(cur_idx) + (m + 1) - d) % (m + 1);
      const T6Slot& org = L.buf[org_idx];
      if (org.timestep == std::numeric_limits<std::int64_t>::min()) continue;

      const std::size_t lag_frames = d * (std::size_t(1) << lev);
      auto it = lag2idx_.find(lag_frames);
      if (it == lag2idx_.end()) continue;
      const std::size_t k = it->second;

      const double dtstep = static_cast<double>(cur.timestep - org.timestep);
      if (dtstep <= 0.0) continue;
      const Tensor6 val = pair_(cur, org);
      update_stats_(k, val, dtstep);
    }

    if (lev + 1 >= levels_.size()) return;

    if (L.pending_idx < 0) {
      L.pending_idx = cur_idx;
      return;
    }

    const T6Slot& a = L.buf[static_cast<std::size_t>(L.pending_idx)];
    const T6Slot& b = L.buf[static_cast<std::size_t>(cur_idx)];
    T6Slot& avg = L.tmp_avg;

    avg_(a, b, avg);
    L.pending_idx = -1;

    Level& N = levels_[lev + 1];
    const std::size_t next_idx = N.write_pos;
    T6Slot& dst = N.buf[next_idx];
    dst.timestep = avg.timestep;

#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
    for (std::size_t p = 0; p < nsel_; ++p) {
      dst.xx[p] = avg.xx[p];
      dst.yy[p] = avg.yy[p];
      dst.zz[p] = avg.zz[p];
      dst.xy[p] = avg.xy[p];
      dst.xz[p] = avg.xz[p];
      dst.yz[p] = avg.yz[p];
    }

    insert_and_accumulate_(lev + 1, static_cast<int>(next_idx));
  }

  static void write_slot_(util::BinaryWriter& w, const T6Slot& s) {
    w.write_i64(s.timestep);
    w.write_u64(static_cast<std::uint64_t>(s.xx.size()));
    if (s.xx.size() != s.yy.size() || s.xx.size() != s.zz.size() || s.xx.size() != s.xy.size() || s.xx.size() != s.xz.size() || s.xx.size() != s.yz.size()) {
      throw std::runtime_error("MultiTauCorrelatorT6: slot channel size mismatch during save_state");
    }
    if (!s.xx.empty()) {
      w.write_bytes(s.xx.data(), sizeof(double) * s.xx.size());
      w.write_bytes(s.yy.data(), sizeof(double) * s.yy.size());
      w.write_bytes(s.zz.data(), sizeof(double) * s.zz.size());
      w.write_bytes(s.xy.data(), sizeof(double) * s.xy.size());
      w.write_bytes(s.xz.data(), sizeof(double) * s.xz.size());
      w.write_bytes(s.yz.data(), sizeof(double) * s.yz.size());
    }
  }

  static void read_slot_(util::BinaryReader& r, T6Slot& s) {
    s.timestep = r.read_i64();
    const std::uint64_t n = r.read_u64();
    if (n > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
      throw std::runtime_error("MultiTauCorrelatorT6: slot size overflow in checkpoint");
    }
    const std::size_t nn = static_cast<std::size_t>(n);
    s.xx.resize(nn);
    s.yy.resize(nn);
    s.zz.resize(nn);
    s.xy.resize(nn);
    s.xz.resize(nn);
    s.yz.resize(nn);
    if (nn > 0) {
      r.read_bytes(s.xx.data(), sizeof(double) * nn);
      r.read_bytes(s.yy.data(), sizeof(double) * nn);
      r.read_bytes(s.zz.data(), sizeof(double) * nn);
      r.read_bytes(s.xy.data(), sizeof(double) * nn);
      r.read_bytes(s.xz.data(), sizeof(double) * nn);
      r.read_bytes(s.yz.data(), sizeof(double) * nn);
    }
  }
};

} // namespace pilots
