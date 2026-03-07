#include <algorithm>
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

using measure_ext::Axis1D;
using measure_ext::axis1d_name;
using measure_ext::axis_length;
using measure_ext::box_volume;
using measure_ext::count_diag_dims;
using measure_ext::dstr;
using measure_ext::get_static_combined_view;
using measure_ext::get_static_group_view;
using measure_ext::lag_axis_name;
using measure_ext::mass_by_atom_from_config;
using measure_ext::orth_area_for_axis;
using measure_ext::parse_axis1d;
using measure_ext::parse_double_list;
using measure_ext::primary_axis_coord;
using measure_ext::resolve_exact_frame_end;
using measure_ext::resolve_path;
using measure_ext::same_index_set;
using measure_ext::shell_volume_3d;
using measure_ext::x_unit_for_axis;

struct RangeOptions {
  std::int64_t frame_start = 0;
  std::int64_t frame_end = -1;
  bool dry_run = false;
};

inline bool frame_in_range(std::size_t frame_index, const RangeOptions& opt) {
  const std::int64_t fi = static_cast<std::int64_t>(frame_index);
  if (fi < opt.frame_start) return false;
  if (opt.frame_end >= 0 && fi > opt.frame_end) return false;
  return true;
}

enum class WeightMode {
  Number,
  Charge,
  Mass,
};

inline WeightMode parse_weight_mode(std::string s) {
  for (auto& c : s) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
  if (s == "number" || s == "count") return WeightMode::Number;
  if (s == "charge" || s == "q") return WeightMode::Charge;
  if (s == "mass") return WeightMode::Mass;
  throw std::runtime_error("invalid weight mode: '" + s + "' (use number|charge|mass)");
}

inline const char* weight_mode_name(WeightMode m) {
  switch (m) {
    case WeightMode::Number: return "number";
    case WeightMode::Charge: return "charge";
    case WeightMode::Mass: return "mass";
  }
  return "number";
}

struct CartesianMSDPairOpT6 {
  Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const {
    const std::size_t n = cur.xx.size();
    if (org.xx.size() != n || cur.yy.size() != n || org.yy.size() != n ||
        cur.zz.size() != n || org.zz.size() != n) {
      throw std::runtime_error("CartesianMSDPairOpT6: slot sizes mismatch");
    }
    double sxx = 0.0, syy = 0.0, szz = 0.0, sxy = 0.0, sxz = 0.0, syz = 0.0;
#if PILOTS_HAS_OPENMP
#pragma omp parallel for reduction(+:sxx,syy,szz,sxy,sxz,syz)
#endif
    for (std::size_t p = 0; p < n; ++p) {
      const double dx = cur.xx[p] - org.xx[p];
      const double dy = cur.yy[p] - org.yy[p];
      const double dz = cur.zz[p] - org.zz[p];
      sxx += dx * dx; syy += dy * dy; szz += dz * dz;
      sxy += dx * dy; sxz += dx * dz; syz += dy * dz;
    }
    const double inv = 1.0 / static_cast<double>(n);
    Tensor6 out;
    out.v[T6_XX] = sxx * inv;
    out.v[T6_YY] = syy * inv;
    out.v[T6_ZZ] = szz * inv;
    out.v[T6_XY] = sxy * inv;
    out.v[T6_XZ] = sxz * inv;
    out.v[T6_YZ] = syz * inv;
    return out;
  }
};

struct DirectionalISFPairOpT6 {
  double q = 1.0;

  Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const {
    const std::size_t n = cur.xx.size();
    if (org.xx.size() != n || cur.yy.size() != n || org.yy.size() != n ||
        cur.zz.size() != n || org.zz.size() != n) {
      throw std::runtime_error("DirectionalISFPairOpT6: slot sizes mismatch");
    }
    double sx = 0.0, sy = 0.0, sz = 0.0;
#if PILOTS_HAS_OPENMP
#pragma omp parallel for reduction(+:sx,sy,sz)
#endif
    for (std::size_t p = 0; p < n; ++p) {
      sx += std::cos(q * (cur.xx[p] - org.xx[p]));
      sy += std::cos(q * (cur.yy[p] - org.yy[p]));
      sz += std::cos(q * (cur.zz[p] - org.zz[p]));
    }
    const double inv = 1.0 / static_cast<double>(n);
    Tensor6 out;
    out.v[T6_XX] = sx * inv;
    out.v[T6_YY] = sy * inv;
    out.v[T6_ZZ] = sz * inv;
    return out;
  }
};

class DirectionalMSDMeasure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    CorrelatorSpec corr;
    bool remove_drift = false;
  };

  DirectionalMSDMeasure(std::string instance_name,
                        std::string output_path,
                        SelectionView sel,
                        SelectionView drift_sel,
                        Options opt)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    drift_name_owned_ = std::string(drift_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    drift_sel_ = SelectionView{drift_name_owned_, drift_sel.idx};
    if (sel_.idx.empty()) throw std::runtime_error("directional_msd: selection is empty");
    corr_ = make_correlator_t6_runtime_auto(sel_.idx.size(), exact_window_frames_(), opt_.corr, CartesianMSDPairOpT6{});
    tmp_xx_.assign(sel_.idx.size(), 0.0);
    tmp_yy_.assign(sel_.idx.size(), 0.0);
    tmp_zz_.assign(sel_.idx.size(), 0.0);
    tmp_zero_.assign(sel_.idx.size(), 0.0);
  }

  std::string type() const override { return "directional_msd"; }
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
    od.columns = {"lag", "time", "msd_xx", "msd_yy", "msd_zz", "msd_xy", "msd_xz", "msd_yz", "msd_trace", "d_xx", "d_yy", "d_zz", "d_trace", "count_pairs", "n_blocks", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    corr_->start();
    started_ = true;
    if (!opt_.range.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");
    double comx = 0.0, comy = 0.0, comz = 0.0;
    if (opt_.remove_drift) {
      double sx = 0.0, sy = 0.0, sz = 0.0;
      for (const std::size_t i : drift_sel_.idx) { sx += xu[i]; sy += yu[i]; sz += zu[i]; }
      const double inv = 1.0 / static_cast<double>(drift_sel_.idx.size());
      comx = sx * inv; comy = sy * inv; comz = sz * inv;
    }
#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
    for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
      const std::size_t i = sel_.idx[p];
      tmp_xx_[p] = xu[i] - comx;
      tmp_yy_[p] = yu[i] - comy;
      tmp_zz_[p] = zu[i] - comz;
    }
    corr_->push(frame.timestep,
                std::span<const double>(tmp_xx_.data(), tmp_xx_.size()),
                std::span<const double>(tmp_yy_.data(), tmp_yy_.size()),
                std::span<const double>(tmp_zz_.data(), tmp_zz_.size()),
                std::span<const double>(tmp_zero_.data(), tmp_zero_.size()),
                std::span<const double>(tmp_zero_.data(), tmp_zero_.size()),
                std::span<const double>(tmp_zero_.data(), tmp_zero_.size()));
  }

  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    const CorrelationSeriesT6 s = corr_->snapshot();
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: directional MSD / Einstein tensor\n";
      ofs << "# columns: lag time msd_xx msd_yy msd_zz msd_xy msd_xz msd_yz msd_trace d_xx d_yy d_zz d_trace count_pairs n_blocks mean_dtimestep\n";
      for (std::size_t i = 0; i < s.lag.size(); ++i) {
        const double trace = s.value_xx[i] + s.value_yy[i] + s.value_zz[i];
        double dxx = 0.0, dyy = 0.0, dzz = 0.0, dtr = 0.0;
        if (s.time[i] > 0.0) {
          dxx = 0.5 * s.value_xx[i] / s.time[i];
          dyy = 0.5 * s.value_yy[i] / s.time[i];
          dzz = 0.5 * s.value_zz[i] / s.time[i];
          dtr = trace / (6.0 * s.time[i]);
        }
        ofs << std::setprecision(17) << s.lag[i] << ' '
            << std::setprecision(17) << s.time[i] << ' '
            << std::setprecision(17) << s.value_xx[i] << ' '
            << std::setprecision(17) << s.value_yy[i] << ' '
            << std::setprecision(17) << s.value_zz[i] << ' '
            << std::setprecision(17) << s.value_xy[i] << ' '
            << std::setprecision(17) << s.value_xz[i] << ' '
            << std::setprecision(17) << s.value_yz[i] << ' '
            << std::setprecision(17) << trace << ' '
            << std::setprecision(17) << dxx << ' '
            << std::setprecision(17) << dyy << ' '
            << std::setprecision(17) << dzz << ' '
            << std::setprecision(17) << dtr << ' '
            << s.count_pairs[i] << ' '
            << s.n_blocks[i] << ' '
            << std::setprecision(17) << s.mean_dtimestep[i] << '\n';
      }
    });
  }

  void finalize() override { if (!opt_.range.dry_run) flush_partial(); }

private:
  std::size_t exact_window_frames_() const {
    if (opt_.corr.type != "exact") return 0;
    if (opt_.range.frame_end < 0) throw std::runtime_error("directional_msd: exact correlator requires finite frame_end");
    return static_cast<std::size_t>(opt_.range.frame_end - opt_.range.frame_start + 1);
  }

  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  std::string drift_name_owned_;
  SelectionView sel_;
  SelectionView drift_sel_;
  Options opt_;
  std::unique_ptr<ICorrelatorT6> corr_;
  bool started_ = false;
  std::vector<double> tmp_xx_, tmp_yy_, tmp_zz_, tmp_zero_;
};

class DirectionalISFMeasure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    CorrelatorSpec corr;
    bool remove_drift = false;
    std::vector<double> q_values;
  };

  DirectionalISFMeasure(std::string instance_name,
                        std::string output_path,
                        SelectionView sel,
                        SelectionView drift_sel,
                        Options opt)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    drift_name_owned_ = std::string(drift_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    drift_sel_ = SelectionView{drift_name_owned_, drift_sel.idx};
    if (sel_.idx.empty()) throw std::runtime_error("directional_isf: selection is empty");
    tmp_xx_.assign(sel_.idx.size(), 0.0);
    tmp_yy_.assign(sel_.idx.size(), 0.0);
    tmp_zz_.assign(sel_.idx.size(), 0.0);
    tmp_zero_.assign(sel_.idx.size(), 0.0);
    for (const double q : opt_.q_values) {
      corr_.push_back(make_correlator_t6_runtime_auto(sel_.idx.size(), exact_window_frames_(), opt_.corr, DirectionalISFPairOpT6{q}));
    }
  }

  std::string type() const override { return "directional_isf"; }
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
    od.columns = {"lag", "time", "q", "fs_x", "fs_y", "fs_z", "fs_mean", "count_pairs", "n_blocks", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));
    md.params["q_values"] = dstr(opt_.q_values.empty() ? 0.0 : opt_.q_values.front());
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    for (auto& c : corr_) c->start();
    started_ = true;
    if (!opt_.range.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");
    double comx = 0.0, comy = 0.0, comz = 0.0;
    if (opt_.remove_drift) {
      double sx = 0.0, sy = 0.0, sz = 0.0;
      for (const std::size_t i : drift_sel_.idx) { sx += xu[i]; sy += yu[i]; sz += zu[i]; }
      const double inv = 1.0 / static_cast<double>(drift_sel_.idx.size());
      comx = sx * inv; comy = sy * inv; comz = sz * inv;
    }
#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
    for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
      const std::size_t i = sel_.idx[p];
      tmp_xx_[p] = xu[i] - comx;
      tmp_yy_[p] = yu[i] - comy;
      tmp_zz_[p] = zu[i] - comz;
    }
    for (auto& c : corr_) {
      c->push(frame.timestep,
              std::span<const double>(tmp_xx_.data(), tmp_xx_.size()),
              std::span<const double>(tmp_yy_.data(), tmp_yy_.size()),
              std::span<const double>(tmp_zz_.data(), tmp_zz_.size()),
              std::span<const double>(tmp_zero_.data(), tmp_zero_.size()),
              std::span<const double>(tmp_zero_.data(), tmp_zero_.size()),
              std::span<const double>(tmp_zero_.data(), tmp_zero_.size()));
    }
  }

  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: directional self-intermediate scattering function\n";
      ofs << "# columns: lag time q fs_x fs_y fs_z fs_mean count_pairs n_blocks mean_dtimestep\n";
      for (std::size_t qi = 0; qi < opt_.q_values.size(); ++qi) {
        const auto s = corr_[qi]->snapshot();
        for (std::size_t i = 0; i < s.lag.size(); ++i) {
          const double mean = (s.value_xx[i] + s.value_yy[i] + s.value_zz[i]) / 3.0;
          ofs << std::setprecision(17) << s.lag[i] << ' '
              << std::setprecision(17) << s.time[i] << ' '
              << std::setprecision(17) << opt_.q_values[qi] << ' '
              << std::setprecision(17) << s.value_xx[i] << ' '
              << std::setprecision(17) << s.value_yy[i] << ' '
              << std::setprecision(17) << s.value_zz[i] << ' '
              << std::setprecision(17) << mean << ' '
              << s.count_pairs[i] << ' '
              << s.n_blocks[i] << ' '
              << std::setprecision(17) << s.mean_dtimestep[i] << '\n';
        }
      }
    });
  }

  void finalize() override { if (!opt_.range.dry_run) flush_partial(); }

private:
  std::size_t exact_window_frames_() const {
    if (opt_.corr.type != "exact") return 0;
    if (opt_.range.frame_end < 0) throw std::runtime_error("directional_isf: exact correlator requires finite frame_end");
    return static_cast<std::size_t>(opt_.range.frame_end - opt_.range.frame_start + 1);
  }

  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  std::string drift_name_owned_;
  SelectionView sel_;
  SelectionView drift_sel_;
  Options opt_;
  bool started_ = false;
  std::vector<std::unique_ptr<ICorrelatorT6>> corr_;
  std::vector<double> tmp_xx_, tmp_yy_, tmp_zz_, tmp_zero_;
};

class SlabRDFMeasure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    Axis1D axis = Axis1D::Z;
    std::size_t n_slab_bins = 0;
    std::size_t n_r_bins = 0;
    double r_max = 0.0;
  };

  SlabRDFMeasure(std::string instance_name,
                 std::string output_path,
                 SelectionView sel_a,
                 SelectionView sel_b,
                 Options opt)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)) {
    sel_a_name_owned_ = std::string(sel_a.name);
    sel_b_name_owned_ = std::string(sel_b.name);
    sel_a_ = SelectionView{sel_a_name_owned_, sel_a.idx};
    sel_b_ = SelectionView{sel_b_name_owned_, sel_b.idx};
    same_sel_ = same_index_set(sel_a_, sel_b_);
    hist_.assign(opt_.n_slab_bins * opt_.n_r_bins, 0.0);
    ref_count_.assign(opt_.n_slab_bins, 0.0);
  }

  std::string type() const override { return "slab_rdf"; }
  std::string instance_name() const override { return instance_name_; }
  output::MeasureDescriptor describe() const override {
    output::MeasureDescriptor md;
    md.instance = instance_name_;
    md.type = type();
    md.selection = std::string(sel_a_.name) + "|" + std::string(sel_b_.name);
    md.n_selected = sel_a_.idx.size() + sel_b_.idx.size();
    output::OutputFileDescriptor od;
    od.path = output_path_;
    od.format = "text";
    od.x_axis = "r";
    od.x_unit = "distance";
    od.columns = {"slab_bin", "slab_center_frac", "r_center", "g", "count_pairs", "n_ref", "n_frames"};
    md.outputs.push_back(std::move(od));
    return md;
  }
  void on_start(const Frame& first_frame) override { (void)first_frame; started_ = true; if (!opt_.range.dry_run) flush_partial(); }
  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");
    const double rho_b = static_cast<double>(sel_b_.idx.size()) / box_volume(frame.box);
    rho_b_sum_ += rho_b;
    for (const std::size_t ia : sel_a_.idx) {
      const double frac = primary_axis_coord(frame.box, xu[ia], yu[ia], zu[ia], opt_.axis) / axis_length(frame.box, opt_.axis);
      std::size_t sb = static_cast<std::size_t>(std::floor(frac * static_cast<double>(opt_.n_slab_bins)));
      if (sb >= opt_.n_slab_bins) sb = opt_.n_slab_bins - 1;
      ref_count_[sb] += 1.0;
      for (const std::size_t jb : sel_b_.idx) {
        if (same_sel_ && ia == jb) continue;
        const auto dr = frame.box.min_image_displacement(xu[ia], yu[ia], zu[ia], xu[jb], yu[jb], zu[jb]);
        const double r = std::sqrt(dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]);
        if (r >= opt_.r_max) continue;
        const std::size_t rb = std::min<std::size_t>(static_cast<std::size_t>(std::floor(r / opt_.r_max * static_cast<double>(opt_.n_r_bins))), opt_.n_r_bins - 1);
        hist_[sb * opt_.n_r_bins + rb] += 1.0;
      }
    }
    ++n_frames_;
  }
  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    const double rho_b = (n_frames_ > 0) ? (rho_b_sum_ / static_cast<double>(n_frames_)) : 0.0;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: slab-resolved RDF using global B density normalization\n";
      ofs << "# axis: " << axis1d_name(opt_.axis) << ", n_slab_bins: " << opt_.n_slab_bins << ", n_r_bins: " << opt_.n_r_bins << "\n";
      ofs << "# columns: slab_bin slab_center_frac r_center g count_pairs n_ref n_frames\n";
      for (std::size_t sb = 0; sb < opt_.n_slab_bins; ++sb) {
        const double slab_center = (static_cast<double>(sb) + 0.5) / static_cast<double>(opt_.n_slab_bins);
        const double nref = ref_count_[sb];
        for (std::size_t rb = 0; rb < opt_.n_r_bins; ++rb) {
          const double r_lo = opt_.r_max * static_cast<double>(rb) / static_cast<double>(opt_.n_r_bins);
          const double r_hi = opt_.r_max * static_cast<double>(rb + 1) / static_cast<double>(opt_.n_r_bins);
          const double shell = shell_volume_3d(r_lo, r_hi);
          const double g = (nref > 0.0 && rho_b > 0.0 && shell > 0.0 && n_frames_ > 0)
              ? (hist_[sb * opt_.n_r_bins + rb] / (nref * rho_b * shell))
              : 0.0;
          ofs << sb << ' ' << std::setprecision(17) << slab_center << ' '
              << std::setprecision(17) << 0.5 * (r_lo + r_hi) << ' '
              << std::setprecision(17) << g << ' '
              << std::setprecision(17) << hist_[sb * opt_.n_r_bins + rb] << ' '
              << std::setprecision(17) << nref << ' '
              << n_frames_ << '\n';
        }
      }
    });
  }
  void finalize() override { if (!opt_.range.dry_run) flush_partial(); }
private:
  std::string instance_name_, output_path_, sel_a_name_owned_, sel_b_name_owned_;
  SelectionView sel_a_, sel_b_;
  Options opt_;
  bool same_sel_ = false;
  bool started_ = false;
  std::vector<double> hist_, ref_count_;
  double rho_b_sum_ = 0.0;
  std::size_t n_frames_ = 0;
};

class CylindricalProfileMeasure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    Axis1D axis = Axis1D::Z;
    std::size_t n_bins = 0;
    double r_max = 0.0;
    WeightMode mode = WeightMode::Number;
    std::string charge_field = "q";
  };

  CylindricalProfileMeasure(std::string instance_name,
                            std::string output_path,
                            SelectionView sel,
                            Options opt,
                            std::vector<double> mass_by_atom)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)),
        mass_by_atom_(std::move(mass_by_atom)) {
    sel_name_owned_ = std::string(sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    weight_sum_.assign(opt_.n_bins, 0.0);
    density_sum_.assign(opt_.n_bins, 0.0);
  }
  std::string type() const override { return "cylindrical_profile"; }
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
    od.x_axis = "r";
    od.x_unit = "distance";
    od.columns = {"bin", "r_center", "mean_weight", "mean_density", "n_frames"};
    md.outputs.push_back(std::move(od));
    return md;
  }
  void on_start(const Frame& first_frame) override { (void)first_frame; started_ = true; if (!opt_.range.dry_run) flush_partial(); }
  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");
    const auto q = (opt_.mode == WeightMode::Charge) ? frame.require_dfield(opt_.charge_field) : std::span<const double>();
    const double axis_len = axis_length(frame.box, opt_.axis);
    const auto center = frame.box.wrap(frame.box.xlo + 0.5 * frame.box.lx(), frame.box.ylo + 0.5 * frame.box.ly(), frame.box.zlo + 0.5 * frame.box.lz());
    for (const std::size_t i : sel_.idx) {
      double w = 1.0;
      if (opt_.mode == WeightMode::Charge) w = q[i];
      else if (opt_.mode == WeightMode::Mass) w = mass_by_atom_.at(i);
      const auto wr = frame.box.wrap(xu[i], yu[i], zu[i]);
      double r = 0.0;
      if (opt_.axis == Axis1D::X) r = std::hypot(wr[1] - center[1], wr[2] - center[2]);
      else if (opt_.axis == Axis1D::Y) r = std::hypot(wr[0] - center[0], wr[2] - center[2]);
      else r = std::hypot(wr[0] - center[0], wr[1] - center[1]);
      if (r >= opt_.r_max) continue;
      const std::size_t b = std::min<std::size_t>(static_cast<std::size_t>(std::floor(r / opt_.r_max * static_cast<double>(opt_.n_bins))), opt_.n_bins - 1);
      const double r_lo = opt_.r_max * static_cast<double>(b) / static_cast<double>(opt_.n_bins);
      const double r_hi = opt_.r_max * static_cast<double>(b + 1) / static_cast<double>(opt_.n_bins);
      const double shell_vol = M_PI * (r_hi * r_hi - r_lo * r_lo) * axis_len;
      weight_sum_[b] += w;
      density_sum_[b] += w / shell_vol;
    }
    ++n_frames_;
  }
  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: cylindrical radial profile\n";
      ofs << "# axis: " << axis1d_name(opt_.axis) << ", mode: " << weight_mode_name(opt_.mode) << "\n";
      ofs << "# columns: bin r_center mean_weight mean_density n_frames\n";
      for (std::size_t b = 0; b < opt_.n_bins; ++b) {
        const double r_lo = opt_.r_max * static_cast<double>(b) / static_cast<double>(opt_.n_bins);
        const double r_hi = opt_.r_max * static_cast<double>(b + 1) / static_cast<double>(opt_.n_bins);
        ofs << b << ' '
            << std::setprecision(17) << 0.5 * (r_lo + r_hi) << ' '
            << std::setprecision(17) << (n_frames_ > 0 ? weight_sum_[b] / static_cast<double>(n_frames_) : 0.0) << ' '
            << std::setprecision(17) << (n_frames_ > 0 ? density_sum_[b] / static_cast<double>(n_frames_) : 0.0) << ' '
            << n_frames_ << '\n';
      }
    });
  }
  void finalize() override { if (!opt_.range.dry_run) flush_partial(); }
private:
  std::string instance_name_, output_path_, sel_name_owned_;
  SelectionView sel_;
  Options opt_;
  std::vector<double> mass_by_atom_;
  bool started_ = false;
  std::vector<double> weight_sum_, density_sum_;
  std::size_t n_frames_ = 0;
};

class Map2DMeasure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    Axis1D axis_u = Axis1D::X;
    Axis1D axis_v = Axis1D::Y;
    std::size_t n_u_bins = 0;
    std::size_t n_v_bins = 0;
    WeightMode mode = WeightMode::Number;
    std::string charge_field = "q";
  };

  Map2DMeasure(std::string instance_name,
               std::string output_path,
               SelectionView sel,
               Options opt,
               std::vector<double> mass_by_atom)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)),
        mass_by_atom_(std::move(mass_by_atom)) {
    sel_name_owned_ = std::string(sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    weight_sum_.assign(opt_.n_u_bins * opt_.n_v_bins, 0.0);
    density_sum_.assign(opt_.n_u_bins * opt_.n_v_bins, 0.0);
  }
  std::string type() const override { return "2d_map"; }
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
    od.x_axis = "u,v";
    od.x_unit = "fractional_coordinate";
    od.columns = {"u_bin", "v_bin", "u_center_frac", "v_center_frac", "mean_weight", "mean_density", "n_frames"};
    md.outputs.push_back(std::move(od));
    return md;
  }
  void on_start(const Frame& first_frame) override { (void)first_frame; started_ = true; if (!opt_.range.dry_run) flush_partial(); }
  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");
    const auto q = (opt_.mode == WeightMode::Charge) ? frame.require_dfield(opt_.charge_field) : std::span<const double>();
    const double vbin = box_volume(frame.box) / static_cast<double>(opt_.n_u_bins * opt_.n_v_bins);
    for (const std::size_t i : sel_.idx) {
      double w = 1.0;
      if (opt_.mode == WeightMode::Charge) w = q[i];
      else if (opt_.mode == WeightMode::Mass) w = mass_by_atom_.at(i);
      const auto lam = frame.box.to_lambda(xu[i], yu[i], zu[i]);
      const auto frac = [&](Axis1D axis) {
        switch (axis) {
          case Axis1D::X: return lam[0] - std::floor(lam[0]);
          case Axis1D::Y: return lam[1] - std::floor(lam[1]);
          case Axis1D::Z: return lam[2] - std::floor(lam[2]);
        }
        return lam[2];
      };
      const std::size_t ub = std::min<std::size_t>(static_cast<std::size_t>(std::floor(frac(opt_.axis_u) * static_cast<double>(opt_.n_u_bins))), opt_.n_u_bins - 1);
      const std::size_t vb = std::min<std::size_t>(static_cast<std::size_t>(std::floor(frac(opt_.axis_v) * static_cast<double>(opt_.n_v_bins))), opt_.n_v_bins - 1);
      weight_sum_[ub * opt_.n_v_bins + vb] += w;
      density_sum_[ub * opt_.n_v_bins + vb] += w / vbin;
    }
    ++n_frames_;
  }
  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: projected 2D map\n";
      ofs << "# axes: " << axis1d_name(opt_.axis_u) << ", " << axis1d_name(opt_.axis_v) << "; mode: " << weight_mode_name(opt_.mode) << "\n";
      ofs << "# columns: u_bin v_bin u_center_frac v_center_frac mean_weight mean_density n_frames\n";
      for (std::size_t ub = 0; ub < opt_.n_u_bins; ++ub) {
        for (std::size_t vb = 0; vb < opt_.n_v_bins; ++vb) {
          const std::size_t idx = ub * opt_.n_v_bins + vb;
          ofs << ub << ' ' << vb << ' '
              << std::setprecision(17) << (static_cast<double>(ub) + 0.5) / static_cast<double>(opt_.n_u_bins) << ' '
              << std::setprecision(17) << (static_cast<double>(vb) + 0.5) / static_cast<double>(opt_.n_v_bins) << ' '
              << std::setprecision(17) << (n_frames_ > 0 ? weight_sum_[idx] / static_cast<double>(n_frames_) : 0.0) << ' '
              << std::setprecision(17) << (n_frames_ > 0 ? density_sum_[idx] / static_cast<double>(n_frames_) : 0.0) << ' '
              << n_frames_ << '\n';
        }
      }
    });
  }
  void finalize() override { if (!opt_.range.dry_run) flush_partial(); }
private:
  std::string instance_name_, output_path_, sel_name_owned_;
  SelectionView sel_;
  Options opt_;
  std::vector<double> mass_by_atom_;
  bool started_ = false;
  std::vector<double> weight_sum_, density_sum_;
  std::size_t n_frames_ = 0;
};

class InterfaceExcessMeasure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    Axis1D axis = Axis1D::Z;
    std::size_t n_bins = 0;
    WeightMode mode = WeightMode::Number;
    std::string charge_field = "q";
    double left_lo_frac = 0.0, left_hi_frac = 0.2;
    double right_lo_frac = 0.8, right_hi_frac = 1.0;
    double divide_frac = 0.5;
  };

  InterfaceExcessMeasure(std::string instance_name,
                         std::string output_path,
                         SelectionView sel,
                         Options opt,
                         std::vector<double> mass_by_atom)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)),
        mass_by_atom_(std::move(mass_by_atom)) {
    sel_name_owned_ = std::string(sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    density_sum_.assign(opt_.n_bins, 0.0);
    axis_length_sum_ = 0.0;
  }
  std::string type() const override { return "interface_excess"; }
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
    od.x_axis = "summary";
    od.x_unit = "none";
    od.columns = {"left_bulk_density", "right_bulk_density", "divide_coord", "gamma_left", "gamma_right", "gamma_total", "n_frames"};
    md.outputs.push_back(std::move(od));
    return md;
  }
  void on_start(const Frame& first_frame) override { (void)first_frame; started_ = true; if (!opt_.range.dry_run) flush_partial(); }
  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");
    const auto q = (opt_.mode == WeightMode::Charge) ? frame.require_dfield(opt_.charge_field) : std::span<const double>();
    const double vbin = box_volume(frame.box) / static_cast<double>(opt_.n_bins);
    const double L = axis_length(frame.box, opt_.axis);
    axis_length_sum_ += L;
    for (const std::size_t i : sel_.idx) {
      double w = 1.0;
      if (opt_.mode == WeightMode::Charge) w = q[i];
      else if (opt_.mode == WeightMode::Mass) w = mass_by_atom_.at(i);
      const double frac = primary_axis_coord(frame.box, xu[i], yu[i], zu[i], opt_.axis) / L;
      std::size_t b = std::min<std::size_t>(static_cast<std::size_t>(std::floor(frac * static_cast<double>(opt_.n_bins))), opt_.n_bins - 1);
      density_sum_[b] += w / vbin;
    }
    ++n_frames_;
  }
  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    const double meanL = (n_frames_ > 0) ? axis_length_sum_ / static_cast<double>(n_frames_) : 0.0;
    auto mean_density = [&](std::size_t b) { return (n_frames_ > 0) ? density_sum_[b] / static_cast<double>(n_frames_) : 0.0; };
    double rhoL = 0.0, rhoR = 0.0;
    std::size_t nL = 0, nR = 0;
    for (std::size_t b = 0; b < opt_.n_bins; ++b) {
      const double c = (static_cast<double>(b) + 0.5) / static_cast<double>(opt_.n_bins);
      if (c >= opt_.left_lo_frac && c <= opt_.left_hi_frac) { rhoL += mean_density(b); ++nL; }
      if (c >= opt_.right_lo_frac && c <= opt_.right_hi_frac) { rhoR += mean_density(b); ++nR; }
    }
    rhoL = (nL > 0) ? rhoL / static_cast<double>(nL) : 0.0;
    rhoR = (nR > 0) ? rhoR / static_cast<double>(nR) : 0.0;
    const double dz = (opt_.n_bins > 0) ? (meanL / static_cast<double>(opt_.n_bins)) : 0.0;
    double gammaL = 0.0, gammaR = 0.0;
    for (std::size_t b = 0; b < opt_.n_bins; ++b) {
      const double c = (static_cast<double>(b) + 0.5) / static_cast<double>(opt_.n_bins);
      if (c < opt_.divide_frac) gammaL += (mean_density(b) - rhoL) * dz;
      else gammaR += (mean_density(b) - rhoR) * dz;
    }
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: interface excess from a 1D profile\n";
      ofs << "# axis: " << axis1d_name(opt_.axis) << ", mode: " << weight_mode_name(opt_.mode) << "\n";
      ofs << "# left_window_frac: [" << opt_.left_lo_frac << ", " << opt_.left_hi_frac << "]\n";
      ofs << "# right_window_frac: [" << opt_.right_lo_frac << ", " << opt_.right_hi_frac << "]\n";
      ofs << "# divide_frac: " << opt_.divide_frac << "\n";
      ofs << "# columns: left_bulk_density right_bulk_density divide_coord gamma_left gamma_right gamma_total n_frames\n";
      ofs << std::setprecision(17) << rhoL << ' '
          << std::setprecision(17) << rhoR << ' '
          << std::setprecision(17) << opt_.divide_frac * meanL << ' '
          << std::setprecision(17) << gammaL << ' '
          << std::setprecision(17) << gammaR << ' '
          << std::setprecision(17) << (gammaL + gammaR) << ' '
          << n_frames_ << '\n';
    });
  }
  void finalize() override { if (!opt_.range.dry_run) flush_partial(); }
private:
  std::string instance_name_, output_path_, sel_name_owned_;
  SelectionView sel_;
  Options opt_;
  std::vector<double> mass_by_atom_;
  bool started_ = false;
  std::vector<double> density_sum_;
  double axis_length_sum_ = 0.0;
  std::size_t n_frames_ = 0;
};

void append_static_caps(const IniConfig& cfg,
                        const std::string& section,
                        MeasureCapabilities& caps) {
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_dfields = {"xu", "yu", "zu"};
  caps.requires_identity_consistent = true;
  caps.scale = ScaleCompatibility{true, true, true};
  caps.group_refs.push_back(cfg.get_string(section, "group", std::optional<std::string>("all")));
  if (cfg.has_key(section, "group_b")) caps.group_refs.push_back(cfg.get_string(section, "group_b"));
  if (cfg.has_key(section, "drift_group")) caps.group_refs.push_back(cfg.get_string(section, "drift_group"));
}

MeasureCapabilities directional_msd_caps(const IniConfig& cfg,
                                         const std::string& section,
                                         const std::string& instance,
                                         const MeasureBuildEnv& env) {
  (void)instance; (void)env; MeasureCapabilities caps; append_static_caps(cfg, section, caps); return caps; }

std::unique_ptr<IMeasure> directional_msd_create(const IniConfig& cfg,
                                                 const std::string& section,
                                                 const std::string& instance,
                                                 const MeasureBuildEnv& env,
                                                 const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("directional_msd factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "directional_msd");
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(false));
  const std::string drift_group = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));
  SelectionView drift_sel = remove_drift ? get_static_group_view(*env.selection_provider, frame0, drift_group, "directional_msd")
                                         : get_static_group_view(*env.selection_provider, frame0, "all", "directional_msd");
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty() ? env.output_dir_general : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("directional_msd.dat"))).lexically_normal();
  DirectionalMSDMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = resolve_exact_frame_end(parse_correlator_spec(cfg, section, env.dt), env.follow,
                                                opt.range.frame_start,
                                                cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1)),
                                                env.input_path);
  opt.range.dry_run = env.dry_run;
  opt.corr = parse_correlator_spec(cfg, section, env.dt);
  opt.remove_drift = remove_drift;
  return std::make_unique<DirectionalMSDMeasure>(instance, out.string(), sel, drift_sel, opt);
}

MeasureCapabilities directional_isf_caps(const IniConfig& cfg,
                                         const std::string& section,
                                         const std::string& instance,
                                         const MeasureBuildEnv& env) {
  (void)instance; (void)env; MeasureCapabilities caps; append_static_caps(cfg, section, caps); return caps; }

std::unique_ptr<IMeasure> directional_isf_create(const IniConfig& cfg,
                                                 const std::string& section,
                                                 const std::string& instance,
                                                 const MeasureBuildEnv& env,
                                                 const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("directional_isf factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "directional_isf");
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(false));
  const std::string drift_group = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));
  SelectionView drift_sel = remove_drift ? get_static_group_view(*env.selection_provider, frame0, drift_group, "directional_isf")
                                         : get_static_group_view(*env.selection_provider, frame0, "all", "directional_isf");
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty() ? env.output_dir_general : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("directional_isf.dat"))).lexically_normal();
  DirectionalISFMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = resolve_exact_frame_end(parse_correlator_spec(cfg, section, env.dt), env.follow,
                                                opt.range.frame_start,
                                                cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1)),
                                                env.input_path);
  opt.range.dry_run = env.dry_run;
  opt.corr = parse_correlator_spec(cfg, section, env.dt);
  opt.remove_drift = remove_drift;
  opt.q_values = parse_double_list(cfg, section, "q_list");
  return std::make_unique<DirectionalISFMeasure>(instance, out.string(), sel, drift_sel, opt);
}

MeasureCapabilities slab_rdf_caps(const IniConfig& cfg,
                                  const std::string& section,
                                  const std::string& instance,
                                  const MeasureBuildEnv& env) {
  (void)instance; (void)env; MeasureCapabilities caps; append_static_caps(cfg, section, caps); return caps; }

std::unique_ptr<IMeasure> slab_rdf_create(const IniConfig& cfg,
                                          const std::string& section,
                                          const std::string& instance,
                                          const MeasureBuildEnv& env,
                                          const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("slab_rdf factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group_a = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo_a = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb_a = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  const std::string group_b = cfg.get_string(section, "group_b", std::optional<std::string>(group_a));
  const std::string topo_b = cfg.get_string(section, "topo_group_b", std::optional<std::string>(topo_a));
  const std::string comb_b = cfg.get_string(section, "combine_b", std::optional<std::string>(comb_a));
  SelectionView sel_a = get_static_combined_view(*env.selection_provider, frame0, group_a, topo_a, comb_a, "slab_rdf");
  SelectionView sel_b = get_static_combined_view(*env.selection_provider, frame0, group_b, topo_b, comb_b, "slab_rdf");
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty() ? env.output_dir_general : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("slab_rdf.dat"))).lexically_normal();
  SlabRDFMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.axis = parse_axis1d(cfg.get_string(section, "axis", std::optional<std::string>("z")));
  opt.n_slab_bins = static_cast<std::size_t>(cfg.get_int64(section, "n_slab_bins", std::optional<std::int64_t>(20)));
  opt.n_r_bins = static_cast<std::size_t>(cfg.get_int64(section, "n_r_bins", std::optional<std::int64_t>(200)));
  opt.r_max = cfg.get_double(section, "r_max");
  return std::make_unique<SlabRDFMeasure>(instance, out.string(), sel_a, sel_b, opt);
}

MeasureCapabilities cylindrical_profile_caps(const IniConfig& cfg,
                                             const std::string& section,
                                             const std::string& instance,
                                             const MeasureBuildEnv& env) {
  (void)instance; (void)env; MeasureCapabilities caps; append_static_caps(cfg, section, caps); const auto mode = parse_weight_mode(cfg.get_string(section, "mode", std::optional<std::string>("number"))); if (mode == WeightMode::Charge) caps.requires_dfields.push_back(cfg.get_string(section, "charge_field", std::optional<std::string>("q"))); if (mode == WeightMode::Mass) { if (cfg.has_key(section, "mass_field")) caps.requires_dfields.push_back(cfg.get_string(section, "mass_field")); else { caps.requires_intfields.push_back("type"); caps.requires_topology_sections.push_back("masses"); } } return caps; }

std::unique_ptr<IMeasure> cylindrical_profile_create(const IniConfig& cfg,
                                                     const std::string& section,
                                                     const std::string& instance,
                                                     const MeasureBuildEnv& env,
                                                     const SystemContext& sysctx) {
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("cylindrical_profile factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "cylindrical_profile");
  const auto mode = parse_weight_mode(cfg.get_string(section, "mode", std::optional<std::string>("number")));
  std::vector<double> mass_by_atom(frame0.natoms, 1.0);
  if (mode == WeightMode::Mass) mass_by_atom = mass_by_atom_from_config(cfg, section, frame0, sysctx);
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty() ? env.output_dir_general : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("cylindrical_profile.dat"))).lexically_normal();
  CylindricalProfileMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.axis = parse_axis1d(cfg.get_string(section, "axis", std::optional<std::string>("z")));
  opt.n_bins = static_cast<std::size_t>(cfg.get_int64(section, "n_bins", std::optional<std::int64_t>(100)));
  opt.r_max = cfg.get_double(section, "r_max");
  opt.mode = mode;
  opt.charge_field = cfg.get_string(section, "charge_field", std::optional<std::string>("q"));
  return std::make_unique<CylindricalProfileMeasure>(instance, out.string(), sel, opt, std::move(mass_by_atom));
}

MeasureCapabilities map2d_caps(const IniConfig& cfg,
                               const std::string& section,
                               const std::string& instance,
                               const MeasureBuildEnv& env) {
  (void)instance; (void)env; MeasureCapabilities caps; append_static_caps(cfg, section, caps); const auto mode = parse_weight_mode(cfg.get_string(section, "mode", std::optional<std::string>("number"))); if (mode == WeightMode::Charge) caps.requires_dfields.push_back(cfg.get_string(section, "charge_field", std::optional<std::string>("q"))); if (mode == WeightMode::Mass) { if (cfg.has_key(section, "mass_field")) caps.requires_dfields.push_back(cfg.get_string(section, "mass_field")); else { caps.requires_intfields.push_back("type"); caps.requires_topology_sections.push_back("masses"); } } return caps; }

std::unique_ptr<IMeasure> map2d_create(const IniConfig& cfg,
                                       const std::string& section,
                                       const std::string& instance,
                                       const MeasureBuildEnv& env,
                                       const SystemContext& sysctx) {
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("2d_map factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "2d_map");
  const auto mode = parse_weight_mode(cfg.get_string(section, "mode", std::optional<std::string>("number")));
  std::vector<double> mass_by_atom(frame0.natoms, 1.0);
  if (mode == WeightMode::Mass) mass_by_atom = mass_by_atom_from_config(cfg, section, frame0, sysctx);
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty() ? env.output_dir_general : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("2d_map.dat"))).lexically_normal();
  Map2DMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.axis_u = parse_axis1d(cfg.get_string(section, "axis_u", std::optional<std::string>("x")));
  opt.axis_v = parse_axis1d(cfg.get_string(section, "axis_v", std::optional<std::string>("y")));
  opt.n_u_bins = static_cast<std::size_t>(cfg.get_int64(section, "n_u_bins", std::optional<std::int64_t>(100)));
  opt.n_v_bins = static_cast<std::size_t>(cfg.get_int64(section, "n_v_bins", std::optional<std::int64_t>(100)));
  opt.mode = mode;
  opt.charge_field = cfg.get_string(section, "charge_field", std::optional<std::string>("q"));
  return std::make_unique<Map2DMeasure>(instance, out.string(), sel, opt, std::move(mass_by_atom));
}

MeasureCapabilities interface_excess_caps(const IniConfig& cfg,
                                          const std::string& section,
                                          const std::string& instance,
                                          const MeasureBuildEnv& env) {
  return map2d_caps(cfg, section, instance, env);
}

std::unique_ptr<IMeasure> interface_excess_create(const IniConfig& cfg,
                                                  const std::string& section,
                                                  const std::string& instance,
                                                  const MeasureBuildEnv& env,
                                                  const SystemContext& sysctx) {
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("interface_excess factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "interface_excess");
  const auto mode = parse_weight_mode(cfg.get_string(section, "mode", std::optional<std::string>("number")));
  std::vector<double> mass_by_atom(frame0.natoms, 1.0);
  if (mode == WeightMode::Mass) mass_by_atom = mass_by_atom_from_config(cfg, section, frame0, sysctx);
  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty() ? env.output_dir_general : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("interface_excess.dat"))).lexically_normal();
  InterfaceExcessMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.axis = parse_axis1d(cfg.get_string(section, "axis", std::optional<std::string>("z")));
  opt.n_bins = static_cast<std::size_t>(cfg.get_int64(section, "n_bins", std::optional<std::int64_t>(200)));
  opt.mode = mode;
  opt.charge_field = cfg.get_string(section, "charge_field", std::optional<std::string>("q"));
  opt.left_lo_frac = cfg.get_double(section, "left_lo_frac", std::optional<double>(0.0));
  opt.left_hi_frac = cfg.get_double(section, "left_hi_frac", std::optional<double>(0.2));
  opt.right_lo_frac = cfg.get_double(section, "right_lo_frac", std::optional<double>(0.8));
  opt.right_hi_frac = cfg.get_double(section, "right_hi_frac", std::optional<double>(1.0));
  opt.divide_frac = cfg.get_double(section, "divide_frac", std::optional<double>(0.5));
  return std::make_unique<InterfaceExcessMeasure>(instance, out.string(), sel, opt, std::move(mass_by_atom));
}

static MeasureRegistrar g_reg_directional_msd("directional_msd", &directional_msd_caps, &directional_msd_create);
static MeasureRegistrar g_reg_directional_isf("directional_isf", &directional_isf_caps, &directional_isf_create);
static MeasureRegistrar g_reg_slab_rdf("slab_rdf", &slab_rdf_caps, &slab_rdf_create);
static MeasureRegistrar g_reg_cylindrical_profile("cylindrical_profile", &cylindrical_profile_caps, &cylindrical_profile_create);
static MeasureRegistrar g_reg_2d_map("2d_map", &map2d_caps, &map2d_create);
static MeasureRegistrar g_reg_interface_excess("interface_excess", &interface_excess_caps, &interface_excess_create);

} // namespace
} // namespace pilots
