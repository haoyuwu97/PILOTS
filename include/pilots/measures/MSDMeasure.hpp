#pragma once

#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <memory>
#include <sstream>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#if PILOTS_HAS_OPENMP
#include <omp.h>
#endif

#include "pilots/correlate/CorrelatorFactory.hpp"
#include "pilots/correlate/CorrelatorSpec.hpp"
#include "pilots/correlate/ICorrelatorT6.hpp"
#include "pilots/correlate/LagAxis.hpp"
#include "pilots/correlate/Tensor6Types.hpp"
#include "pilots/core/Frame.hpp"
#include "pilots/measures/IMeasure.hpp"
#include "pilots/select/Selection.hpp"
#include "pilots/select/SelectionView.hpp"
#include "pilots/util/AtomicFile.hpp"
#include "pilots/util/BinaryIO.hpp"

namespace pilots {

// MSD pair operator on a vector-valued signal stored in (xx,yy,zz) channels.
// - Input mapping: x := slot.xx, y := slot.yy, z := slot.zz
// - Output: a scalar MSD value stored into the xx component of Tensor6
//   (other channels are 0). This preserves correct SEM for the chosen combination.
//
// diag_mask bits: xx=1, yy=2, zz=4.
struct MSDPairOpT6 {
  int diag_mask = 7;

  Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const {
    const bool use_x = (diag_mask & 1) != 0;
    const bool use_y = (diag_mask & 2) != 0;
    const bool use_z = (diag_mask & 4) != 0;

    const std::size_t n = cur.xx.size();
    if (org.xx.size() != n || cur.yy.size() != n || org.yy.size() != n || cur.zz.size() != n || org.zz.size() != n) {
      throw std::runtime_error("MSDPairOpT6: slot sizes mismatch");
    }

    double sum_dr2 = 0.0;
#if PILOTS_HAS_OPENMP
#pragma omp parallel for reduction(+:sum_dr2)
#endif
    for (std::size_t p = 0; p < n; ++p) {
      const double dx = cur.xx[p] - org.xx[p];
      const double dy = cur.yy[p] - org.yy[p];
      const double dz = cur.zz[p] - org.zz[p];
      double dr2 = 0.0;
      if (use_x) dr2 += dx * dx;
      if (use_y) dr2 += dy * dy;
      if (use_z) dr2 += dz * dz;
      sum_dr2 += dr2;
    }

    Tensor6 out;
    out.v[T6_XX] = sum_dr2 / static_cast<double>(n);
    return out;
  }
};

class MSDMeasure final : public IMeasure {
public:
  struct Options {
    std::int64_t frame_start = 0;
    std::int64_t frame_end = -1; // inclusive; for correlator=exact must be finite

    int diag_mask = 7;      // xx=1,yy=2,zz=4
    bool remove_drift = true;

    CorrelatorSpec corr;    // shared correlator configuration
  };

  MSDMeasure(std::string instance_name, std::string output_path, SelectionView sel, SelectionView drift_sel, Options opt)
  : instance_name_(std::move(instance_name)), output_path_(std::move(output_path)), opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    drift_name_owned_ = std::string(drift_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    drift_sel_ = SelectionView{drift_name_owned_, drift_sel.idx};
    if (opt_.frame_start < 0) throw std::runtime_error("MSDMeasure: frame_start must be >= 0");
    if (opt_.frame_end >= 0 && opt_.frame_end < opt_.frame_start) {
      throw std::runtime_error("MSDMeasure: frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.corr.dt > 0.0)) throw std::runtime_error("MSDMeasure: corr.dt must be > 0");
    if (opt_.diag_mask <= 0 || opt_.diag_mask > 7) throw std::runtime_error("MSDMeasure: invalid diag_mask");
    if (sel_.idx.empty()) throw std::runtime_error("MSDMeasure: selection is empty");
    if (opt_.remove_drift && drift_sel_.idx.empty()) throw std::runtime_error("MSDMeasure: drift selection is empty");
    if (opt_.corr.axis == LagAxis::TimeBin && !(opt_.corr.timebin_width > 0.0)) {
      throw std::runtime_error("MSDMeasure: timebin_width must be > 0 for lag_axis=timebin");
    }

    if (opt_.corr.type == "exact") {
      if (opt_.frame_end < 0) throw std::runtime_error("MSDMeasure: correlator=exact requires finite frame_end");
    }

    std::size_t window_frames = 0;
    if (opt_.corr.type == "exact") {
      window_frames = static_cast<std::size_t>(opt_.frame_end - opt_.frame_start + 1);
      if (window_frames < 2) throw std::runtime_error("MSDMeasure: exact window must have >=2 frames");
    }

    // Measures only provide a PairOp object; the factory performs type-erasure
    // and supplies a default multi-tau AvgOp (element-wise) when needed.
    const MSDPairOpT6 pair_op{opt_.diag_mask};
    corr_ = make_correlator_t6_runtime_auto(sel_.idx.size(), window_frames, opt_.corr, pair_op);
  }

  std::string type() const override { return "msd"; }
  std::string instance_name() const override { return instance_name_; }

  output::MeasureDescriptor describe() const override {
    output::MeasureDescriptor md;
    md.instance = instance_name_;
    md.type = type();
    md.selection = std::string(sel_.name);
    md.n_selected = sel_.idx.size();

    output::OutputFileDescriptor od;
    od.path = output_path_;
    od.format = "text";
    // Dataset axis semantics for downstream tooling (results.json).
    od.x_axis = lag_axis_name(opt_.corr.axis);
    if (opt_.corr.axis == LagAxis::Frame) {
      od.x_unit = "frames";
    } else if (opt_.corr.axis == LagAxis::Timestep) {
      od.x_unit = "timesteps";
    } else if (opt_.corr.axis == LagAxis::TimeBin) {
      od.x_unit = "simulation_time";
    }
    od.columns = {"lag", "time", "msd", "count_pairs", "sem", "n_blocks", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));

    auto dstr = [](double v) {
      std::ostringstream oss;
      oss << std::setprecision(17) << v;
      return oss.str();
    };
    md.params["correlator"] = opt_.corr.type;
    md.params["lag_axis"] = lag_axis_name(opt_.corr.axis);
    md.params["dt"] = dstr(opt_.corr.dt);
    md.params["frame_start"] = std::to_string(opt_.frame_start);
    md.params["frame_end"] = std::to_string(opt_.frame_end);
    md.params["diag_mask"] = std::to_string(opt_.diag_mask);
    md.params["remove_drift"] = opt_.remove_drift ? "1" : "0";
    md.params["drift_group"] = drift_sel_.name;
    md.params["block_size"] = std::to_string(opt_.corr.block_size);
    if (opt_.corr.type == "exact") {
      md.params["lag_stride"] = std::to_string(opt_.corr.lag_stride);
    }
    if (opt_.corr.type == "multitau") {
      md.params["mt_channels"] = std::to_string(opt_.corr.mt_channels);
      md.params["mt_levels"] = std::to_string(opt_.corr.mt_levels);
    }
    if (opt_.corr.axis == LagAxis::TimeBin) {
      md.params["timebin_width"] = dstr(opt_.corr.timebin_width);
      md.params["min_pairs_per_bin"] = std::to_string(opt_.corr.min_pairs_per_bin);
      md.params["bin_merge"] = opt_.corr.bin_merge ? "1" : "0";
    }
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;

    const std::size_t nsel = sel_.idx.size();
    tmp_xx_.assign(nsel, 0.0);
    tmp_yy_.assign(nsel, 0.0);
    tmp_zz_.assign(nsel, 0.0);
    tmp_xy_.assign(nsel, 0.0);
    tmp_xz_.assign(nsel, 0.0);
    tmp_yz_.assign(nsel, 0.0);

    corr_->start();
    started_ = true;

    // Create an initial output file (header-only) so follow/online mode has something to tail.
    flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    const std::int64_t fi = static_cast<std::int64_t>(frame_index);
    if (fi < opt_.frame_start) return;
    if (opt_.frame_end >= 0 && fi > opt_.frame_end) return;

    // Strictly require core coordinate fields.
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");

    // Drift COM over drift_sel
    double comx = 0.0, comy = 0.0, comz = 0.0;
    if (opt_.remove_drift) {
      const std::size_t nd = drift_sel_.idx.size();
      double sx = 0.0, sy = 0.0, sz = 0.0;
#if PILOTS_HAS_OPENMP
#pragma omp parallel for reduction(+:sx,sy,sz)
#endif
      for (std::size_t k = 0; k < nd; ++k) {
        const std::size_t i = drift_sel_.idx[k];
        sx += xu[i];
        sy += yu[i];
        sz += zu[i];
      }
      const double inv = 1.0 / static_cast<double>(nd);
      comx = sx * inv;
      comy = sy * inv;
      comz = sz * inv;
    }

    const std::size_t nsel = sel_.idx.size();
#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
    for (std::size_t p = 0; p < nsel; ++p) {
      const std::size_t i = sel_.idx[p];
      tmp_xx_[p] = xu[i] - comx;
      tmp_yy_[p] = yu[i] - comy;
      tmp_zz_[p] = zu[i] - comz;
    }
    // tmp_xy/tmp_xz/tmp_yz remain zero

    corr_->push(frame.timestep,
               std::span<const double>(tmp_xx_.data(), tmp_xx_.size()),
               std::span<const double>(tmp_yy_.data(), tmp_yy_.size()),
               std::span<const double>(tmp_zz_.data(), tmp_zz_.size()),
               std::span<const double>(tmp_xy_.data(), tmp_xy_.size()),
               std::span<const double>(tmp_xz_.data(), tmp_xz_.size()),
               std::span<const double>(tmp_yz_.data(), tmp_yz_.size()));
  }

  void flush_partial() override {
    if (!started_) return;
    const CorrelationSeriesT6 series = corr_->snapshot();

    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      write_header_(ofs);
      write_series_(ofs, series);
    });
  }

  void finalize() override {
    // finalize is non-destructive, but represents the final flush.
    flush_partial();
  }

  void save_state(std::ostream& os) const override {
    util::write_magic(os, "PILOTSMSDMEASURE");
    util::BinaryWriter w(os);
    w.write_u32(1);

    w.write_i64(opt_.frame_start);
    w.write_i64(opt_.frame_end);
    w.write_i32(opt_.diag_mask);
    w.write_u8(opt_.remove_drift ? 1 : 0);

    w.write_string(opt_.corr.type);
    w.write_i32(static_cast<std::int32_t>(opt_.corr.axis));
    w.write_f64(opt_.corr.dt);
    w.write_f64(opt_.corr.timebin_width);
    w.write_u64(opt_.corr.min_pairs_per_bin);
    w.write_u8(opt_.corr.bin_merge ? 1 : 0);
    w.write_u64(static_cast<std::uint64_t>(opt_.corr.lag_stride));
    w.write_u64(static_cast<std::uint64_t>(opt_.corr.block_size));
    w.write_u64(static_cast<std::uint64_t>(opt_.corr.mt_channels));
    w.write_u64(static_cast<std::uint64_t>(opt_.corr.mt_levels));

    w.write_u64(static_cast<std::uint64_t>(sel_.idx.size()));
    w.write_u64(static_cast<std::uint64_t>(drift_sel_.idx.size()));
    w.write_u8(started_ ? 1 : 0);

    // Correlator state (may be large)
    corr_->save_state(os);

    util::write_magic(os, "PILOTSEND");
  }

  void load_state(std::istream& is) override {
    util::require_magic(is, "PILOTSMSDMEASURE");
    util::BinaryReader r(is);
    const std::uint32_t ver = r.read_u32();
    if (ver != 1) throw std::runtime_error("MSDMeasure: unsupported checkpoint version");

    auto near_eq = [](double a, double b) {
      const double d = (a > b) ? (a - b) : (b - a);
      return d <= 1e-12;
    };

    const std::int64_t fs = r.read_i64();
    const std::int64_t fe = r.read_i64();
    const std::int32_t dm = r.read_i32();
    const bool rd = (r.read_u8() != 0);

    if (fs != opt_.frame_start || fe != opt_.frame_end) throw std::runtime_error("MSDMeasure: frame range mismatch in checkpoint");
    if (dm != opt_.diag_mask) throw std::runtime_error("MSDMeasure: diag_mask mismatch in checkpoint");
    if (rd != opt_.remove_drift) throw std::runtime_error("MSDMeasure: remove_drift mismatch in checkpoint");

    const std::string ct = r.read_string();
    const LagAxis ax = static_cast<LagAxis>(r.read_i32());
    const double dt = r.read_f64();
    const double tb = r.read_f64();
    const std::uint64_t mp = r.read_u64();
    const bool bm = (r.read_u8() != 0);
    const std::uint64_t ls = r.read_u64();
    const std::uint64_t bs = r.read_u64();
    const std::uint64_t mch = r.read_u64();
    const std::uint64_t mlv = r.read_u64();

    if (ct != opt_.corr.type) throw std::runtime_error("MSDMeasure: correlator type mismatch in checkpoint");
    if (ax != opt_.corr.axis) throw std::runtime_error("MSDMeasure: lag_axis mismatch in checkpoint");
    if (!near_eq(dt, opt_.corr.dt)) throw std::runtime_error("MSDMeasure: dt mismatch in checkpoint");
    if (!near_eq(tb, opt_.corr.timebin_width)) throw std::runtime_error("MSDMeasure: timebin_width mismatch in checkpoint");
    if (mp != opt_.corr.min_pairs_per_bin) throw std::runtime_error("MSDMeasure: min_pairs_per_bin mismatch in checkpoint");
    if (bm != opt_.corr.bin_merge) throw std::runtime_error("MSDMeasure: bin_merge mismatch in checkpoint");
    if (ls != static_cast<std::uint64_t>(opt_.corr.lag_stride)) throw std::runtime_error("MSDMeasure: lag_stride mismatch in checkpoint");
    if (bs != static_cast<std::uint64_t>(opt_.corr.block_size)) throw std::runtime_error("MSDMeasure: block_size mismatch in checkpoint");
    if (mch != static_cast<std::uint64_t>(opt_.corr.mt_channels)) throw std::runtime_error("MSDMeasure: mt_channels mismatch in checkpoint");
    if (mlv != static_cast<std::uint64_t>(opt_.corr.mt_levels)) throw std::runtime_error("MSDMeasure: mt_levels mismatch in checkpoint");

    const std::uint64_t nsel = r.read_u64();
    const std::uint64_t nd = r.read_u64();
    if (nsel != static_cast<std::uint64_t>(sel_.idx.size())) throw std::runtime_error("MSDMeasure: selection size mismatch in checkpoint");
    if (nd != static_cast<std::uint64_t>(drift_sel_.idx.size())) throw std::runtime_error("MSDMeasure: drift selection size mismatch in checkpoint");

    started_ = (r.read_u8() != 0);

    corr_->load_state(is);
    util::require_magic(is, "PILOTSEND");

    // After loading state, ensure the on-disk output is consistent.
    flush_partial();
  }

private:
  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  std::string drift_name_owned_;
  SelectionView sel_;
  SelectionView drift_sel_;
  Options opt_;

  std::unique_ptr<ICorrelatorT6> corr_;
  bool started_ = false;

  // scratch buffers
  std::vector<double> tmp_xx_, tmp_yy_, tmp_zz_, tmp_xy_, tmp_xz_, tmp_yz_;

  void write_header_(std::ostream& ofs) const {
    ofs << "# PILOTS: MSD (tensor6 correlator; scalar stored in xx channel)\n";
    ofs << "# correlator: " << opt_.corr.type << "\n";
    ofs << "# lag_axis: " << lag_axis_name(opt_.corr.axis) << "\n";
    if (opt_.corr.axis == LagAxis::TimeBin) {
      ofs << "# timebin_width: " << std::setprecision(17) << opt_.corr.timebin_width << "\n";
      ofs << "# min_pairs_per_bin: " << opt_.corr.min_pairs_per_bin << "\n";
      ofs << "# bin_merge: " << (opt_.corr.bin_merge ? 1 : 0) << "\n";
    }
    ofs << "# dt: " << std::setprecision(17) << opt_.corr.dt << " (seconds per timestep)\n";
    ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
    ofs << "# remove_drift: " << (opt_.remove_drift ? 1 : 0) << " (drift_group=" << drift_sel_.name << ", n=" << drift_sel_.idx.size() << ")\n";
    ofs << "# frame_range: [" << opt_.frame_start << ", " << opt_.frame_end << "]\n";
    ofs << "# components_mask: " << opt_.diag_mask << " (xx=1,yy=2,zz=4)\n";
    if (opt_.corr.type == "exact") {
      ofs << "# exact_lag_stride: " << opt_.corr.lag_stride << "\n";
    } else if (opt_.corr.type == "multitau") {
      ofs << "# multitau_channels: " << opt_.corr.mt_channels << ", levels: " << opt_.corr.mt_levels << "\n";
    }
    ofs << "# block_size: " << opt_.corr.block_size << " (0 disables SEM)\n";
    ofs << "# columns: lag  time  msd  count_pairs  sem  n_blocks  mean_dtimestep\n";
  }

  static void write_series_(std::ostream& ofs, const CorrelationSeriesT6& series) {
    const std::size_t M = series.lag.size();
    for (std::size_t i = 0; i < M; ++i) {
      ofs << std::setprecision(17) << series.lag[i] << " "
          << std::setprecision(17) << series.time[i] << " "
          << std::setprecision(17) << series.value_xx[i] << " "
          << series.count_pairs[i] << " "
          << std::setprecision(17) << series.sem_xx[i] << " "
          << series.n_blocks[i] << " "
          << std::setprecision(17) << series.mean_dtimestep[i]
          << "\n";
    }
  }
};

} // namespace pilots
