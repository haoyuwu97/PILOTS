#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#if PILOTS_HAS_OPENMP
#include <omp.h>
#endif

#include "MeasureCommon.hpp"
#include "pilots/correlate/CorrelatorFactory.hpp"
#include "pilots/correlate/ICorrelatorT6.hpp"
#include "pilots/correlate/Tensor6Types.hpp"
#include "pilots/measures/IMeasure.hpp"
#include "pilots/measures/MeasureRegistry.hpp"
#include "pilots/select/SelectionView.hpp"
#include "pilots/util/AtomicFile.hpp"

namespace fs = std::filesystem;

namespace pilots {
namespace {

using measure_ext::SelectedChains;
using measure_ext::append_integer_like_field_cap;
using measure_ext::OrderedChains;
using measure_ext::chain_id_per_atom_from_config;
using measure_ext::build_selected_chains;
using measure_ext::count_diag_dims;
using measure_ext::dstr;
using measure_ext::get_static_combined_view;
using measure_ext::get_static_group_view;
using measure_ext::join_ints;
using measure_ext::ordered_chains_from_config;
using measure_ext::parse_diag_mask;
using measure_ext::parse_int_list;
using measure_ext::resolve_exact_frame_end;
using measure_ext::resolve_measure_output_path;
using measure_ext::x_unit_for_axis;

constexpr double kPi = 3.141592653589793238462643383279502884;

enum class GMode {
  G1,
  G2,
  G3,
};

inline const char* gmode_name(GMode m) {
  switch (m) {
    case GMode::G1: return "g1";
    case GMode::G2: return "g2";
    case GMode::G3: return "g3";
  }
  return "g1";
}

struct PolymerDynamicsOptions {
  std::int64_t frame_start = 0;
  std::int64_t frame_end = -1;
  int diag_mask = 7;
  bool remove_drift = true;
  bool dry_run = false;
  CorrelatorSpec corr;
  GMode mode = GMode::G1;
};

struct SimpleMSDPairOpT6 {
  int diag_mask = 7;

  Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const {
    const bool use_x = (diag_mask & 1) != 0;
    const bool use_y = (diag_mask & 2) != 0;
    const bool use_z = (diag_mask & 4) != 0;

    const std::size_t n = cur.xx.size();
    if (org.xx.size() != n || cur.yy.size() != n || org.yy.size() != n ||
        cur.zz.size() != n || org.zz.size() != n) {
      throw std::runtime_error("SimpleMSDPairOpT6: slot sizes mismatch");
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

struct ModeAutoCorrPairOpT6 {
  int diag_mask = 7;

  Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const {
    const bool use_x = (diag_mask & 1) != 0;
    const bool use_y = (diag_mask & 2) != 0;
    const bool use_z = (diag_mask & 4) != 0;
    const std::size_t n = cur.xx.size();
    if (org.xx.size() != n || cur.yy.size() != n || org.yy.size() != n ||
        cur.zz.size() != n || org.zz.size() != n) {
      throw std::runtime_error("ModeAutoCorrPairOpT6: slot sizes mismatch");
    }

    double sum = 0.0;
#if PILOTS_HAS_OPENMP
#pragma omp parallel for reduction(+:sum)
#endif
    for (std::size_t p = 0; p < n; ++p) {
      double dot = 0.0;
      if (use_x) dot += cur.xx[p] * org.xx[p];
      if (use_y) dot += cur.yy[p] * org.yy[p];
      if (use_z) dot += cur.zz[p] * org.zz[p];
      sum += dot;
    }

    Tensor6 out;
    out.v[T6_XX] = sum / static_cast<double>(n);
    return out;
  }
};

class PolymerGMeasure final : public IMeasure {
public:
  PolymerGMeasure(std::string instance_name,
                  std::string output_path,
                  SelectionView sel,
                  SelectionView drift_sel,
                  PolymerDynamicsOptions opt,
                  SelectedChains chains)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)),
        chains_(std::move(chains)) {
    sel_name_owned_ = std::string(sel.name);
    drift_name_owned_ = std::string(drift_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    drift_sel_ = SelectionView{drift_name_owned_, drift_sel.idx};

    validate_common_();
    if (opt_.mode != GMode::G1 && chains_.n_chains() == 0) {
      throw std::runtime_error(std::string(gmode_name(opt_.mode)) + ": no chains are available in the selected atoms");
    }

    inv_chain_size_.resize(chains_.n_chains(), 0.0);
    for (std::size_t c = 0; c < chains_.n_chains(); ++c) {
      if (chains_.selpos_by_chain[c].empty()) {
        throw std::runtime_error(std::string(gmode_name(opt_.mode)) + ": an empty chain was constructed");
      }
      inv_chain_size_[c] = 1.0 / static_cast<double>(chains_.selpos_by_chain[c].size());
    }

    const std::size_t nsignal = (opt_.mode == GMode::G3) ? chains_.n_chains() : sel_.idx.size();
    corr_ = make_correlator_t6_runtime_auto(nsignal, exact_window_frames_(), opt_.corr, SimpleMSDPairOpT6{opt_.diag_mask});
  }

  std::string type() const override { return gmode_name(opt_.mode); }
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
    od.x_axis = lag_axis_name(opt_.corr.axis);
    od.x_unit = x_unit_for_axis(opt_.corr.axis);
    od.columns = {"lag", "time", type(), "count_pairs", "sem", "n_blocks", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));

    fill_common_params_(md);
    md.params["n_chains"] = std::to_string(chains_.n_chains());
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    tmp_sel_x_.assign(sel_.idx.size(), 0.0);
    tmp_sel_y_.assign(sel_.idx.size(), 0.0);
    tmp_sel_z_.assign(sel_.idx.size(), 0.0);
    tmp_chain_comx_.assign(chains_.n_chains(), 0.0);
    tmp_chain_comy_.assign(chains_.n_chains(), 0.0);
    tmp_chain_comz_.assign(chains_.n_chains(), 0.0);

    const std::size_t nsignal = (opt_.mode == GMode::G3) ? chains_.n_chains() : sel_.idx.size();
    tmp_out_x_.assign(nsignal, 0.0);
    tmp_out_y_.assign(nsignal, 0.0);
    tmp_out_z_.assign(nsignal, 0.0);
    tmp_out_xy_.assign(nsignal, 0.0);
    tmp_out_xz_.assign(nsignal, 0.0);
    tmp_out_yz_.assign(nsignal, 0.0);

    corr_->start();
    started_ = true;
    if (!opt_.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range_(frame_index)) return;
    fill_selected_positions_(frame);

    if (opt_.mode == GMode::G1) {
      corr_->push(frame.timestep,
                  std::span<const double>(tmp_sel_x_.data(), tmp_sel_x_.size()),
                  std::span<const double>(tmp_sel_y_.data(), tmp_sel_y_.size()),
                  std::span<const double>(tmp_sel_z_.data(), tmp_sel_z_.size()),
                  std::span<const double>(tmp_out_xy_.data(), tmp_out_xy_.size()),
                  std::span<const double>(tmp_out_xz_.data(), tmp_out_xz_.size()),
                  std::span<const double>(tmp_out_yz_.data(), tmp_out_yz_.size()));
      return;
    }

    compute_chain_coms_();

    if (opt_.mode == GMode::G2) {
#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
      for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
        const std::size_t c = chains_.chain_of_selpos[p];
        tmp_out_x_[p] = tmp_sel_x_[p] - tmp_chain_comx_[c];
        tmp_out_y_[p] = tmp_sel_y_[p] - tmp_chain_comy_[c];
        tmp_out_z_[p] = tmp_sel_z_[p] - tmp_chain_comz_[c];
      }
      corr_->push(frame.timestep,
                  std::span<const double>(tmp_out_x_.data(), tmp_out_x_.size()),
                  std::span<const double>(tmp_out_y_.data(), tmp_out_y_.size()),
                  std::span<const double>(tmp_out_z_.data(), tmp_out_z_.size()),
                  std::span<const double>(tmp_out_xy_.data(), tmp_out_xy_.size()),
                  std::span<const double>(tmp_out_xz_.data(), tmp_out_xz_.size()),
                  std::span<const double>(tmp_out_yz_.data(), tmp_out_yz_.size()));
      return;
    }

    for (std::size_t c = 0; c < chains_.n_chains(); ++c) {
      tmp_out_x_[c] = tmp_chain_comx_[c];
      tmp_out_y_[c] = tmp_chain_comy_[c];
      tmp_out_z_[c] = tmp_chain_comz_[c];
    }
    corr_->push(frame.timestep,
                std::span<const double>(tmp_out_x_.data(), tmp_out_x_.size()),
                std::span<const double>(tmp_out_y_.data(), tmp_out_y_.size()),
                std::span<const double>(tmp_out_z_.data(), tmp_out_z_.size()),
                std::span<const double>(tmp_out_xy_.data(), tmp_out_xy_.size()),
                std::span<const double>(tmp_out_xz_.data(), tmp_out_xz_.size()),
                std::span<const double>(tmp_out_yz_.data(), tmp_out_yz_.size()));
  }

  void flush_partial() override {
    if (opt_.dry_run || !started_) return;
    const CorrelationSeriesT6 series = corr_->snapshot();
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      write_header_(ofs);
      for (std::size_t i = 0; i < series.lag.size(); ++i) {
        ofs << std::setprecision(17) << series.lag[i] << ' '
            << std::setprecision(17) << series.time[i] << ' '
            << std::setprecision(17) << series.value_xx[i] << ' '
            << series.count_pairs[i] << ' '
            << std::setprecision(17) << series.sem_xx[i] << ' '
            << series.n_blocks[i] << ' '
            << std::setprecision(17) << series.mean_dtimestep[i] << '\n';
      }
    });
  }

  void finalize() override {
    if (opt_.dry_run) return;
    flush_partial();
  }

private:
  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  std::string drift_name_owned_;
  SelectionView sel_;
  SelectionView drift_sel_;
  PolymerDynamicsOptions opt_;
  SelectedChains chains_;
  std::unique_ptr<ICorrelatorT6> corr_;
  bool started_ = false;

  std::vector<double> inv_chain_size_;
  std::vector<double> tmp_sel_x_, tmp_sel_y_, tmp_sel_z_;
  std::vector<double> tmp_chain_comx_, tmp_chain_comy_, tmp_chain_comz_;
  std::vector<double> tmp_out_x_, tmp_out_y_, tmp_out_z_, tmp_out_xy_, tmp_out_xz_, tmp_out_yz_;

  void validate_common_() const {
    if (opt_.frame_start < 0) throw std::runtime_error(std::string(gmode_name(opt_.mode)) + ": frame_start must be >= 0");
    if (opt_.frame_end >= 0 && opt_.frame_end < opt_.frame_start) {
      throw std::runtime_error(std::string(gmode_name(opt_.mode)) + ": frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.corr.dt > 0.0)) throw std::runtime_error(std::string(gmode_name(opt_.mode)) + ": corr.dt must be > 0");
    if (opt_.diag_mask <= 0 || opt_.diag_mask > 7) throw std::runtime_error(std::string(gmode_name(opt_.mode)) + ": invalid diag_mask");
    if (sel_.idx.empty()) throw std::runtime_error(std::string(gmode_name(opt_.mode)) + ": selection is empty");
    if (opt_.remove_drift && drift_sel_.idx.empty()) throw std::runtime_error(std::string(gmode_name(opt_.mode)) + ": drift selection is empty");
    if (opt_.corr.axis == LagAxis::TimeBin && !(opt_.corr.timebin_width > 0.0)) {
      throw std::runtime_error(std::string(gmode_name(opt_.mode)) + ": timebin_width must be > 0 for lag_axis=timebin");
    }
    if (opt_.corr.type == "exact" && opt_.frame_end < 0) {
      throw std::runtime_error(std::string(gmode_name(opt_.mode)) + ": correlator=exact requires finite frame_end");
    }
  }

  std::size_t exact_window_frames_() const {
    if (opt_.corr.type != "exact") return 0;
    const std::size_t w = static_cast<std::size_t>(opt_.frame_end - opt_.frame_start + 1);
    if (w < 2) throw std::runtime_error(std::string(gmode_name(opt_.mode)) + ": exact window must have >= 2 frames");
    return w;
  }

  bool frame_in_range_(std::size_t frame_index) const {
    const std::int64_t fi = static_cast<std::int64_t>(frame_index);
    if (fi < opt_.frame_start) return false;
    if (opt_.frame_end >= 0 && fi > opt_.frame_end) return false;
    return true;
  }

  void fill_common_params_(output::MeasureDescriptor& md) const {
    md.params["correlator"] = opt_.corr.type;
    md.params["lag_axis"] = lag_axis_name(opt_.corr.axis);
    md.params["dt"] = dstr(opt_.corr.dt);
    md.params["frame_start"] = std::to_string(opt_.frame_start);
    md.params["frame_end"] = std::to_string(opt_.frame_end);
    md.params["diag_mask"] = std::to_string(opt_.diag_mask);
    md.params["dimensions"] = std::to_string(count_diag_dims(opt_.diag_mask));
    md.params["remove_drift"] = opt_.remove_drift ? "1" : "0";
    md.params["drift_group"] = std::string(drift_sel_.name);
    md.params["block_size"] = std::to_string(opt_.corr.block_size);
    if (opt_.corr.type == "exact") md.params["lag_stride"] = std::to_string(opt_.corr.lag_stride);
    if (opt_.corr.type == "multitau") {
      md.params["mt_channels"] = std::to_string(opt_.corr.mt_channels);
      md.params["mt_levels"] = std::to_string(opt_.corr.mt_levels);
    }
    if (opt_.corr.axis == LagAxis::TimeBin) {
      md.params["timebin_width"] = dstr(opt_.corr.timebin_width);
      md.params["min_pairs_per_bin"] = std::to_string(opt_.corr.min_pairs_per_bin);
      md.params["bin_merge"] = opt_.corr.bin_merge ? "1" : "0";
    }
  }

  void write_header_(std::ostream& ofs) const {
    ofs << "# PILOTS: polymer " << type() << "\n";
    ofs << "# correlator: " << opt_.corr.type << "\n";
    ofs << "# lag_axis: " << lag_axis_name(opt_.corr.axis) << "\n";
    if (opt_.corr.axis == LagAxis::TimeBin) {
      ofs << "# timebin_width: " << std::setprecision(17) << opt_.corr.timebin_width << "\n";
      ofs << "# min_pairs_per_bin: " << opt_.corr.min_pairs_per_bin << "\n";
      ofs << "# bin_merge: " << (opt_.corr.bin_merge ? 1 : 0) << "\n";
    }
    ofs << "# dt: " << std::setprecision(17) << opt_.corr.dt << " (seconds per timestep)\n";
    ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
    ofs << "# n_chains: " << chains_.n_chains() << "\n";
    ofs << "# remove_drift: " << (opt_.remove_drift ? 1 : 0)
        << " (drift_group=" << drift_sel_.name << ", n=" << drift_sel_.idx.size() << ")\n";
    ofs << "# frame_range: [" << opt_.frame_start << ", " << opt_.frame_end << "]\n";
    ofs << "# components_mask: " << opt_.diag_mask << " (xx=1,yy=2,zz=4)\n";
    ofs << "# dimensions: " << count_diag_dims(opt_.diag_mask) << "\n";
    if (opt_.corr.type == "exact") {
      ofs << "# exact_lag_stride: " << opt_.corr.lag_stride << "\n";
    } else if (opt_.corr.type == "multitau") {
      ofs << "# multitau_channels: " << opt_.corr.mt_channels << ", levels: " << opt_.corr.mt_levels << "\n";
    }
    ofs << "# block_size: " << opt_.corr.block_size << " (0 disables SEM)\n";
    ofs << "# columns: lag  time  " << type() << "  count_pairs  sem  n_blocks  mean_dtimestep\n";
  }

  void fill_selected_positions_(const Frame& frame) {
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");

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

#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
    for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
      const std::size_t i = sel_.idx[p];
      tmp_sel_x_[p] = xu[i] - comx;
      tmp_sel_y_[p] = yu[i] - comy;
      tmp_sel_z_[p] = zu[i] - comz;
    }
  }

  void compute_chain_coms_() {
    for (std::size_t c = 0; c < chains_.n_chains(); ++c) {
      double sx = 0.0, sy = 0.0, sz = 0.0;
      for (const std::size_t p : chains_.selpos_by_chain[c]) {
        sx += tmp_sel_x_[p];
        sy += tmp_sel_y_[p];
        sz += tmp_sel_z_[p];
      }
      tmp_chain_comx_[c] = sx * inv_chain_size_[c];
      tmp_chain_comy_[c] = sy * inv_chain_size_[c];
      tmp_chain_comz_[c] = sz * inv_chain_size_[c];
    }
  }
};

class RouseMeasure final : public IMeasure {
public:
  struct ModeSpec {
    int p = 1;
    std::vector<std::size_t> chain_slots;
    std::vector<std::vector<double>> weights;
    std::unique_ptr<ICorrelatorT6> corr;
    std::vector<double> x, y, z, xy, xz, yz;
  };

  RouseMeasure(std::string instance_name,
               std::string output_path,
               SelectionView sel,
               SelectionView drift_sel,
               PolymerDynamicsOptions opt,
               OrderedChains chains,
               std::vector<int> p_values)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)),
        chains_(std::move(chains)),
        p_values_(std::move(p_values)) {
    sel_name_owned_ = std::string(sel.name);
    drift_name_owned_ = std::string(drift_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    drift_sel_ = SelectionView{drift_name_owned_, drift_sel.idx};

    validate_common_();
    if (chains_.n_chains() == 0) {
      throw std::runtime_error("rouse: no ordered chains are available in the selected atoms");
    }

    build_modes_();
  }

  std::string type() const override { return "rouse"; }
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
    od.x_axis = lag_axis_name(opt_.corr.axis);
    od.x_unit = x_unit_for_axis(opt_.corr.axis);
    od.columns = {"lag", "time", "p", "cp", "cp_norm", "count_pairs", "sem_cp", "n_blocks", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));

    fill_common_params_(md);
    md.params["p_values"] = join_ints(p_values_);
    md.params["n_ordered_chains"] = std::to_string(chains_.n_chains());
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    tmp_sel_x_.assign(sel_.idx.size(), 0.0);
    tmp_sel_y_.assign(sel_.idx.size(), 0.0);
    tmp_sel_z_.assign(sel_.idx.size(), 0.0);
    for (auto& m : modes_) {
      m.x.assign(m.chain_slots.size(), 0.0);
      m.y.assign(m.chain_slots.size(), 0.0);
      m.z.assign(m.chain_slots.size(), 0.0);
      m.xy.assign(m.chain_slots.size(), 0.0);
      m.xz.assign(m.chain_slots.size(), 0.0);
      m.yz.assign(m.chain_slots.size(), 0.0);
      m.corr->start();
    }
    started_ = true;
    if (!opt_.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range_(frame_index)) return;
    fill_selected_positions_(frame);

    for (auto& m : modes_) {
      for (std::size_t ic = 0; ic < m.chain_slots.size(); ++ic) {
        const std::size_t c = m.chain_slots[ic];
        const auto& chain = chains_.selpos_by_chain[c];
        const auto& w = m.weights[ic];
        double sx = 0.0, sy = 0.0, sz = 0.0;
        for (std::size_t n = 0; n < chain.size(); ++n) {
          const std::size_t p = chain[n];
          sx += w[n] * tmp_sel_x_[p];
          sy += w[n] * tmp_sel_y_[p];
          sz += w[n] * tmp_sel_z_[p];
        }
        m.x[ic] = sx;
        m.y[ic] = sy;
        m.z[ic] = sz;
      }

      m.corr->push(frame.timestep,
                   std::span<const double>(m.x.data(), m.x.size()),
                   std::span<const double>(m.y.data(), m.y.size()),
                   std::span<const double>(m.z.data(), m.z.size()),
                   std::span<const double>(m.xy.data(), m.xy.size()),
                   std::span<const double>(m.xz.data(), m.xz.size()),
                   std::span<const double>(m.yz.data(), m.yz.size()));
    }
  }

  void flush_partial() override {
    if (opt_.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      write_header_(ofs);
      for (const auto& m : modes_) {
        const CorrelationSeriesT6 series = m.corr->snapshot();
        const double c0 = (!series.value_xx.empty()) ? series.value_xx.front() : 0.0;
        for (std::size_t i = 0; i < series.lag.size(); ++i) {
          const double cp = series.value_xx[i];
          const double cpn = (std::abs(c0) > 1e-300) ? (cp / c0) : 0.0;
          ofs << std::setprecision(17) << series.lag[i] << ' '
              << std::setprecision(17) << series.time[i] << ' '
              << m.p << ' '
              << std::setprecision(17) << cp << ' '
              << std::setprecision(17) << cpn << ' '
              << series.count_pairs[i] << ' '
              << std::setprecision(17) << series.sem_xx[i] << ' '
              << series.n_blocks[i] << ' '
              << std::setprecision(17) << series.mean_dtimestep[i] << '\n';
        }
      }
    });
  }

  void finalize() override {
    if (opt_.dry_run) return;
    flush_partial();
  }

private:
  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  std::string drift_name_owned_;
  SelectionView sel_;
  SelectionView drift_sel_;
  PolymerDynamicsOptions opt_;
  OrderedChains chains_;
  std::vector<int> p_values_;
  std::vector<ModeSpec> modes_;
  bool started_ = false;

  std::vector<double> tmp_sel_x_, tmp_sel_y_, tmp_sel_z_;

  void validate_common_() const {
    if (opt_.frame_start < 0) throw std::runtime_error("rouse: frame_start must be >= 0");
    if (opt_.frame_end >= 0 && opt_.frame_end < opt_.frame_start) {
      throw std::runtime_error("rouse: frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.corr.dt > 0.0)) throw std::runtime_error("rouse: corr.dt must be > 0");
    if (opt_.diag_mask <= 0 || opt_.diag_mask > 7) throw std::runtime_error("rouse: invalid diag_mask");
    if (sel_.idx.empty()) throw std::runtime_error("rouse: selection is empty");
    if (opt_.remove_drift && drift_sel_.idx.empty()) throw std::runtime_error("rouse: drift selection is empty");
    if (opt_.corr.axis == LagAxis::TimeBin && !(opt_.corr.timebin_width > 0.0)) {
      throw std::runtime_error("rouse: timebin_width must be > 0 for lag_axis=timebin");
    }
    if (opt_.corr.type == "exact" && opt_.frame_end < 0) {
      throw std::runtime_error("rouse: correlator=exact requires finite frame_end");
    }
  }

  std::size_t exact_window_frames_() const {
    if (opt_.corr.type != "exact") return 0;
    const std::size_t w = static_cast<std::size_t>(opt_.frame_end - opt_.frame_start + 1);
    if (w < 2) throw std::runtime_error("rouse: exact window must have >= 2 frames");
    return w;
  }

  bool frame_in_range_(std::size_t frame_index) const {
    const std::int64_t fi = static_cast<std::int64_t>(frame_index);
    if (fi < opt_.frame_start) return false;
    if (opt_.frame_end >= 0 && fi > opt_.frame_end) return false;
    return true;
  }

  void build_modes_() {
    const std::size_t window = exact_window_frames_();
    modes_.clear();
    modes_.reserve(p_values_.size());

    for (const int p : p_values_) {
      if (p < 0) throw std::runtime_error("rouse: p_values must be >= 0");

      ModeSpec m;
      m.p = p;
      for (std::size_t c = 0; c < chains_.n_chains(); ++c) {
        const std::size_t n = chains_.selpos_by_chain[c].size();
        if (n == 0) continue;
        if (p > 0 && n <= static_cast<std::size_t>(p)) continue;

        m.chain_slots.push_back(c);
        m.weights.emplace_back();
        auto& w = m.weights.back();
        w.resize(n, 0.0);
        const double invn = 1.0 / static_cast<double>(n);
        for (std::size_t k = 0; k < n; ++k) {
          w[k] = invn * std::cos(kPi * static_cast<double>(p) * (static_cast<double>(k) + 0.5) / static_cast<double>(n));
        }
      }

      if (m.chain_slots.empty()) {
        throw std::runtime_error("rouse: requested p=" + std::to_string(p) + " has no valid chains in the selection");
      }

      m.corr = make_correlator_t6_runtime_auto(
          m.chain_slots.size(), window, opt_.corr, ModeAutoCorrPairOpT6{opt_.diag_mask});
      modes_.push_back(std::move(m));
    }
  }

  void fill_selected_positions_(const Frame& frame) {
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");

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

#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
    for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
      const std::size_t i = sel_.idx[p];
      tmp_sel_x_[p] = xu[i] - comx;
      tmp_sel_y_[p] = yu[i] - comy;
      tmp_sel_z_[p] = zu[i] - comz;
    }
  }

  void fill_common_params_(output::MeasureDescriptor& md) const {
    md.params["correlator"] = opt_.corr.type;
    md.params["lag_axis"] = lag_axis_name(opt_.corr.axis);
    md.params["dt"] = dstr(opt_.corr.dt);
    md.params["frame_start"] = std::to_string(opt_.frame_start);
    md.params["frame_end"] = std::to_string(opt_.frame_end);
    md.params["diag_mask"] = std::to_string(opt_.diag_mask);
    md.params["dimensions"] = std::to_string(count_diag_dims(opt_.diag_mask));
    md.params["remove_drift"] = opt_.remove_drift ? "1" : "0";
    md.params["drift_group"] = std::string(drift_sel_.name);
    md.params["block_size"] = std::to_string(opt_.corr.block_size);
    if (opt_.corr.type == "exact") md.params["lag_stride"] = std::to_string(opt_.corr.lag_stride);
    if (opt_.corr.type == "multitau") {
      md.params["mt_channels"] = std::to_string(opt_.corr.mt_channels);
      md.params["mt_levels"] = std::to_string(opt_.corr.mt_levels);
    }
    if (opt_.corr.axis == LagAxis::TimeBin) {
      md.params["timebin_width"] = dstr(opt_.corr.timebin_width);
      md.params["min_pairs_per_bin"] = std::to_string(opt_.corr.min_pairs_per_bin);
      md.params["bin_merge"] = opt_.corr.bin_merge ? "1" : "0";
    }
  }

  void write_header_(std::ostream& ofs) const {
    ofs << "# PILOTS: Rouse-mode autocorrelation\n";
    ofs << "# correlator: " << opt_.corr.type << "\n";
    ofs << "# lag_axis: " << lag_axis_name(opt_.corr.axis) << "\n";
    if (opt_.corr.axis == LagAxis::TimeBin) {
      ofs << "# timebin_width: " << std::setprecision(17) << opt_.corr.timebin_width << "\n";
      ofs << "# min_pairs_per_bin: " << opt_.corr.min_pairs_per_bin << "\n";
      ofs << "# bin_merge: " << (opt_.corr.bin_merge ? 1 : 0) << "\n";
    }
    ofs << "# dt: " << std::setprecision(17) << opt_.corr.dt << " (seconds per timestep)\n";
    ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
    ofs << "# n_ordered_chains: " << chains_.n_chains() << "\n";
    ofs << "# p_values: " << join_ints(p_values_) << "\n";
    ofs << "# remove_drift: " << (opt_.remove_drift ? 1 : 0)
        << " (drift_group=" << drift_sel_.name << ", n=" << drift_sel_.idx.size() << ")\n";
    ofs << "# frame_range: [" << opt_.frame_start << ", " << opt_.frame_end << "]\n";
    ofs << "# components_mask: " << opt_.diag_mask << " (xx=1,yy=2,zz=4)\n";
    ofs << "# dimensions: " << count_diag_dims(opt_.diag_mask) << "\n";
    if (opt_.corr.type == "exact") {
      ofs << "# exact_lag_stride: " << opt_.corr.lag_stride << "\n";
    } else if (opt_.corr.type == "multitau") {
      ofs << "# multitau_channels: " << opt_.corr.mt_channels << ", levels: " << opt_.corr.mt_levels << "\n";
    }
    for (const auto& m : modes_) {
      ofs << "# active_chains[p=" << m.p << "]: " << m.chain_slots.size() << "\n";
    }
    ofs << "# block_size: " << opt_.corr.block_size << " (0 disables SEM)\n";
    ofs << "# columns: lag  time  p  cp  cp_norm  count_pairs  sem_cp  n_blocks  mean_dtimestep\n";
  }
};

void append_common_dynamics_caps(const IniConfig& cfg,
                                 const std::string& section,
                                 MeasureCapabilities& caps) {
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_identity_consistent = true;
  caps.requires_dfields = {"xu", "yu", "zu"};
  caps.scale = ScaleCompatibility{true, true, true};
  const std::string group_ref = cfg.get_string(section, "group", std::optional<std::string>("all"));
  caps.group_refs.push_back(group_ref);
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(true));
  if (remove_drift) {
    const std::string drift_group_ref = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));
    caps.group_refs.push_back(drift_group_ref);
  }
}

void append_chain_membership_caps(const IniConfig& cfg,
                                  const std::string& section,
                                  const MeasureBuildEnv& env,
                                  MeasureCapabilities& caps) {
  if (cfg.has_key(section, "chain_id_field")) {
    append_integer_like_field_cap(caps, cfg.get_string(section, "chain_id_field"));
    return;
  }
  if (env.first_frame && env.first_frame->has_mol) {
    caps.requires_i64fields.push_back("mol");
    return;
  }
  caps.requires_topology_sections.push_back("bonds");
}

void append_chain_order_caps(const IniConfig& cfg,
                             const std::string& section,
                             const MeasureBuildEnv& env,
                             MeasureCapabilities& caps) {
  append_chain_membership_caps(cfg, section, env, caps);
  if (cfg.has_key(section, "chain_pos_field")) {
    caps.requires_dfields.push_back(cfg.get_string(section, "chain_pos_field"));
  } else {
    caps.requires_topology_sections.push_back("bonds");
  }
}

PolymerDynamicsOptions make_polymer_dynamics_options(const IniConfig& cfg,
                                                     const std::string& section,
                                                     const MeasureBuildEnv& env,
                                                     GMode mode) {
  const CorrelatorSpec corr_spec = parse_correlator_spec(cfg, section, env.dt);
  const std::int64_t frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  const std::int64_t frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  if (frame_start < 0) throw std::runtime_error("frame_start must be >= 0");
  if (frame_end >= 0 && frame_end < frame_start) {
    throw std::runtime_error("frame_end must be -1 or >= frame_start");
  }
  const int diag_mask = parse_diag_mask(cfg.get_string(section, "components", std::optional<std::string>("xxyyzz")));

  PolymerDynamicsOptions opt;
  opt.frame_start = frame_start;
  opt.frame_end = (corr_spec.type == "exact")
      ? resolve_exact_frame_end(corr_spec, env.follow, frame_start, frame_end, env.input_path)
      : frame_end;
  opt.diag_mask = diag_mask;
  opt.remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(true));
  opt.dry_run = env.dry_run;
  opt.corr = corr_spec;
  opt.mode = mode;
  return opt;
}

std::pair<SelectionView, SelectionView> resolve_polymer_dynamics_selections(const IniConfig& cfg,
                                                                            const std::string& section,
                                                                            const MeasureBuildEnv& env,
                                                                            const char* measure_name) {
  if (!env.first_frame) throw std::runtime_error(std::string(measure_name) + " factory: first_frame is null");
  if (!env.selection_provider) throw std::runtime_error(std::string(measure_name) + " factory: SelectionProvider is missing");
  const Frame& frame0 = *env.first_frame;
  const std::string group_ref = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo_group_ref = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string combine_expr = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(true));
  const std::string drift_group_ref = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));

  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group_ref, topo_group_ref, combine_expr, measure_name);
  SelectionView drift_sel = remove_drift
      ? get_static_group_view(*env.selection_provider, frame0, drift_group_ref, std::string(measure_name) + " drift_group")
      : get_static_group_view(*env.selection_provider, frame0, "all", std::string(measure_name) + " drift_group");
  return {sel, drift_sel};
}

fs::path resolve_polymer_output_path(const IniConfig& cfg,
                                     const std::string& section,
                                     const MeasureBuildEnv& env,
                                     const std::string& default_name) {
  return resolve_measure_output_path(cfg, section, env, "output", default_name);
}

MeasureCapabilities g1_caps(const IniConfig& cfg,
                            const std::string& section,
                            const std::string& instance,
                            const MeasureBuildEnv& env) {
  (void)instance;
  MeasureCapabilities caps;
  append_common_dynamics_caps(cfg, section, caps);
  (void)env;
  return caps;
}

std::unique_ptr<IMeasure> g1_create(const IniConfig& cfg,
                                    const std::string& section,
                                    const std::string& instance,
                                    const MeasureBuildEnv& env,
                                    const SystemContext& sysctx) {
  (void)sysctx;
  auto [sel, drift_sel] = resolve_polymer_dynamics_selections(cfg, section, env, "g1");
  const PolymerDynamicsOptions opt = make_polymer_dynamics_options(cfg, section, env, GMode::G1);
  const fs::path out_path = resolve_polymer_output_path(cfg, section, env, "g1.dat");
  return std::make_unique<PolymerGMeasure>(instance, out_path.string(), sel, drift_sel, opt, SelectedChains{});
}

MeasureCapabilities g2_caps(const IniConfig& cfg,
                            const std::string& section,
                            const std::string& instance,
                            const MeasureBuildEnv& env) {
  (void)instance;
  MeasureCapabilities caps;
  append_common_dynamics_caps(cfg, section, caps);
  append_chain_membership_caps(cfg, section, env, caps);
  return caps;
}

std::unique_ptr<IMeasure> g2_create(const IniConfig& cfg,
                                    const std::string& section,
                                    const std::string& instance,
                                    const MeasureBuildEnv& env,
                                    const SystemContext& sysctx) {
  auto [sel, drift_sel] = resolve_polymer_dynamics_selections(cfg, section, env, "g2");
  if (!env.first_frame) throw std::runtime_error("g2 factory: first_frame is null");
  const auto chain_id = chain_id_per_atom_from_config(cfg, section, *env.first_frame, sysctx);
  const SelectedChains chains = build_selected_chains(sel, std::span<const std::int64_t>(chain_id.data(), chain_id.size()));
  const PolymerDynamicsOptions opt = make_polymer_dynamics_options(cfg, section, env, GMode::G2);
  const fs::path out_path = resolve_polymer_output_path(cfg, section, env, "g2.dat");
  return std::make_unique<PolymerGMeasure>(instance, out_path.string(), sel, drift_sel, opt, chains);
}

MeasureCapabilities g3_caps(const IniConfig& cfg,
                            const std::string& section,
                            const std::string& instance,
                            const MeasureBuildEnv& env) {
  (void)instance;
  MeasureCapabilities caps;
  append_common_dynamics_caps(cfg, section, caps);
  append_chain_membership_caps(cfg, section, env, caps);
  return caps;
}

std::unique_ptr<IMeasure> g3_create(const IniConfig& cfg,
                                    const std::string& section,
                                    const std::string& instance,
                                    const MeasureBuildEnv& env,
                                    const SystemContext& sysctx) {
  auto [sel, drift_sel] = resolve_polymer_dynamics_selections(cfg, section, env, "g3");
  if (!env.first_frame) throw std::runtime_error("g3 factory: first_frame is null");
  const auto chain_id = chain_id_per_atom_from_config(cfg, section, *env.first_frame, sysctx);
  const SelectedChains chains = build_selected_chains(sel, std::span<const std::int64_t>(chain_id.data(), chain_id.size()));
  const PolymerDynamicsOptions opt = make_polymer_dynamics_options(cfg, section, env, GMode::G3);
  const fs::path out_path = resolve_polymer_output_path(cfg, section, env, "g3.dat");
  return std::make_unique<PolymerGMeasure>(instance, out_path.string(), sel, drift_sel, opt, chains);
}

MeasureCapabilities rouse_caps(const IniConfig& cfg,
                               const std::string& section,
                               const std::string& instance,
                               const MeasureBuildEnv& env) {
  (void)instance;
  MeasureCapabilities caps;
  append_common_dynamics_caps(cfg, section, caps);
  append_chain_order_caps(cfg, section, env, caps);
  return caps;
}

std::vector<int> parse_p_values(const IniConfig& cfg, const std::string& section) {
  if (cfg.has_key(section, "p_values")) {
    return parse_int_list(cfg, section, "p_values");
  }
  const int pmax = static_cast<int>(cfg.get_int64(section, "p_max", std::optional<std::int64_t>(4)));
  if (pmax < 0) throw std::runtime_error("rouse: p_max must be >= 0");
  std::vector<int> out;
  out.reserve(static_cast<std::size_t>(pmax) + 1);
  for (int p = 1; p <= pmax; ++p) out.push_back(p);
  if (out.empty()) out.push_back(0);
  return out;
}

std::unique_ptr<IMeasure> rouse_create(const IniConfig& cfg,
                                       const std::string& section,
                                       const std::string& instance,
                                       const MeasureBuildEnv& env,
                                       const SystemContext& sysctx) {
  auto [sel, drift_sel] = resolve_polymer_dynamics_selections(cfg, section, env, "rouse");
  if (!env.first_frame) throw std::runtime_error("rouse factory: first_frame is null");
  const auto chain_id = chain_id_per_atom_from_config(cfg, section, *env.first_frame, sysctx);
  const SelectedChains groups = build_selected_chains(sel, std::span<const std::int64_t>(chain_id.data(), chain_id.size()));
  const OrderedChains chains = ordered_chains_from_config(cfg, section, *env.first_frame, sysctx, sel, groups);
  const PolymerDynamicsOptions opt = make_polymer_dynamics_options(cfg, section, env, GMode::G1);
  const std::vector<int> p_values = parse_p_values(cfg, section);
  const fs::path out_path = resolve_polymer_output_path(cfg, section, env, "rouse.dat");
  return std::make_unique<RouseMeasure>(instance, out_path.string(), sel, drift_sel, opt, chains, p_values);
}

static MeasureRegistrar g_register_g1("g1", &g1_caps, &g1_create);
static MeasureRegistrar g_register_g2("g2", &g2_caps, &g2_create);
static MeasureRegistrar g_register_g3("g3", &g3_caps, &g3_create);
static MeasureRegistrar g_register_rouse("rouse", &rouse_caps, &rouse_create);

} // namespace
} // namespace pilots
