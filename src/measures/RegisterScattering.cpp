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

using measure_ext::count_diag_dims;
using measure_ext::dstr;
using measure_ext::get_static_combined_view;
using measure_ext::get_static_group_view;
using measure_ext::join_doubles;
using measure_ext::parse_diag_mask;
using measure_ext::parse_double_list;
using measure_ext::resolve_exact_frame_end;
using measure_ext::resolve_path;
using measure_ext::x_unit_for_axis;

inline double spherical_bessel_j0(double x) {
  if (std::abs(x) < 1e-14) return 1.0;
  return std::sin(x) / x;
}

struct SelfISFPairOpT6 {
  double q = 1.0;
  int diag_mask = 7;

  Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const {
    const bool use_x = (diag_mask & 1) != 0;
    const bool use_y = (diag_mask & 2) != 0;
    const bool use_z = (diag_mask & 4) != 0;
    const std::size_t n = cur.xx.size();
    if (org.xx.size() != n || cur.yy.size() != n || org.yy.size() != n ||
        cur.zz.size() != n || org.zz.size() != n) {
      throw std::runtime_error("SelfISFPairOpT6: slot sizes mismatch");
    }

    double sum = 0.0;
#if PILOTS_HAS_OPENMP
#pragma omp parallel for reduction(+:sum)
#endif
    for (std::size_t p = 0; p < n; ++p) {
      const double dx = cur.xx[p] - org.xx[p];
      const double dy = cur.yy[p] - org.yy[p];
      const double dz = cur.zz[p] - org.zz[p];
      double dr2 = 0.0;
      if (use_x) dr2 += dx * dx;
      if (use_y) dr2 += dy * dy;
      if (use_z) dr2 += dz * dz;
      sum += spherical_bessel_j0(q * std::sqrt(dr2));
    }

    Tensor6 out;
    out.v[T6_XX] = sum / static_cast<double>(n);
    return out;
  }
};

struct Chi4OverlapPairOpT6 {
  double a = 1.0;
  int diag_mask = 7;

  Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const {
    const bool use_x = (diag_mask & 1) != 0;
    const bool use_y = (diag_mask & 2) != 0;
    const bool use_z = (diag_mask & 4) != 0;
    const double a2 = a * a;
    const std::size_t n = cur.xx.size();
    if (org.xx.size() != n || cur.yy.size() != n || org.yy.size() != n ||
        cur.zz.size() != n || org.zz.size() != n) {
      throw std::runtime_error("Chi4OverlapPairOpT6: slot sizes mismatch");
    }

    double qsum = 0.0;
#if PILOTS_HAS_OPENMP
#pragma omp parallel for reduction(+:qsum)
#endif
    for (std::size_t p = 0; p < n; ++p) {
      const double dx = cur.xx[p] - org.xx[p];
      const double dy = cur.yy[p] - org.yy[p];
      const double dz = cur.zz[p] - org.zz[p];
      double dr2 = 0.0;
      if (use_x) dr2 += dx * dx;
      if (use_y) dr2 += dy * dy;
      if (use_z) dr2 += dz * dz;
      qsum += (dr2 <= a2) ? 1.0 : 0.0;
    }

    const double qmean = qsum / static_cast<double>(n);
    Tensor6 out;
    out.v[T6_XX] = qmean;
    out.v[T6_YY] = qmean * qmean;
    return out;
  }
};

struct Chi4ISFPairOpT6 {
  double q = 1.0;
  int diag_mask = 7;

  Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const {
    const bool use_x = (diag_mask & 1) != 0;
    const bool use_y = (diag_mask & 2) != 0;
    const bool use_z = (diag_mask & 4) != 0;
    const std::size_t n = cur.xx.size();
    if (org.xx.size() != n || cur.yy.size() != n || org.yy.size() != n ||
        cur.zz.size() != n || org.zz.size() != n) {
      throw std::runtime_error("Chi4ISFPairOpT6: slot sizes mismatch");
    }

    double fsum = 0.0;
#if PILOTS_HAS_OPENMP
#pragma omp parallel for reduction(+:fsum)
#endif
    for (std::size_t p = 0; p < n; ++p) {
      const double dx = cur.xx[p] - org.xx[p];
      const double dy = cur.yy[p] - org.yy[p];
      const double dz = cur.zz[p] - org.zz[p];
      double dr2 = 0.0;
      if (use_x) dr2 += dx * dx;
      if (use_y) dr2 += dy * dy;
      if (use_z) dr2 += dz * dz;
      fsum += spherical_bessel_j0(q * std::sqrt(dr2));
    }

    const double fmean = fsum / static_cast<double>(n);
    Tensor6 out;
    out.v[T6_XX] = fmean;
    out.v[T6_YY] = fmean * fmean;
    return out;
  }
};

struct FamilyCommonOptions {
  std::int64_t frame_start = 0;
  std::int64_t frame_end = -1;
  int diag_mask = 7;
  bool remove_drift = true;
  bool dry_run = false;
  CorrelatorSpec corr;
};

class SelfISFMeasure final : public IMeasure {
public:
  SelfISFMeasure(std::string instance_name,
                 std::string output_path,
                 SelectionView sel,
                 SelectionView drift_sel,
                 FamilyCommonOptions opt,
                 std::vector<double> q_values)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)),
        q_values_(std::move(q_values)) {
    sel_name_owned_ = std::string(sel.name);
    drift_name_owned_ = std::string(drift_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    drift_sel_ = SelectionView{drift_name_owned_, drift_sel.idx};

    validate_common_("SelfISFMeasure");
    for (const double q : q_values_) {
      if (!(q > 0.0)) throw std::runtime_error("SelfISFMeasure: each q must be > 0");
      corr_.push_back(make_correlator_t6_runtime_auto(
          sel_.idx.size(), exact_window_frames_(), opt_.corr, SelfISFPairOpT6{q, opt_.diag_mask}));
    }
  }

  std::string type() const override { return "self_isf"; }
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
    od.columns = {"lag", "time", "q", "fs", "count_pairs", "sem", "n_blocks", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));

    fill_common_params_(md);
    md.params["q_values"] = join_doubles(q_values_);
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    allocate_scratch_();
    for (auto& c : corr_) c->start();
    started_ = true;
    if (!opt_.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range_(frame_index)) return;
    fill_centered_positions_(frame);
    for (auto& c : corr_) {
      c->push(frame.timestep,
              std::span<const double>(tmp_xx_.data(), tmp_xx_.size()),
              std::span<const double>(tmp_yy_.data(), tmp_yy_.size()),
              std::span<const double>(tmp_zz_.data(), tmp_zz_.size()),
              std::span<const double>(tmp_xy_.data(), tmp_xy_.size()),
              std::span<const double>(tmp_xz_.data(), tmp_xz_.size()),
              std::span<const double>(tmp_yz_.data(), tmp_yz_.size()));
    }
  }

  void flush_partial() override {
    if (opt_.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      write_common_header_(ofs, "self intermediate scattering function", "q_values", join_doubles(q_values_));
      ofs << "# columns: lag  time  q  fs  count_pairs  sem  n_blocks  mean_dtimestep\n";
      for (std::size_t iq = 0; iq < q_values_.size(); ++iq) {
        const CorrelationSeriesT6 series = corr_[iq]->snapshot();
        for (std::size_t i = 0; i < series.lag.size(); ++i) {
          ofs << std::setprecision(17) << series.lag[i] << ' '
              << std::setprecision(17) << series.time[i] << ' '
              << std::setprecision(17) << q_values_[iq] << ' '
              << std::setprecision(17) << series.value_xx[i] << ' '
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
  FamilyCommonOptions opt_;
  std::vector<double> q_values_;
  std::vector<std::unique_ptr<ICorrelatorT6>> corr_;
  bool started_ = false;

  std::vector<double> tmp_xx_, tmp_yy_, tmp_zz_, tmp_xy_, tmp_xz_, tmp_yz_;

  void validate_common_(const char* tag) const {
    if (opt_.frame_start < 0) throw std::runtime_error(std::string(tag) + ": frame_start must be >= 0");
    if (opt_.frame_end >= 0 && opt_.frame_end < opt_.frame_start) {
      throw std::runtime_error(std::string(tag) + ": frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.corr.dt > 0.0)) throw std::runtime_error(std::string(tag) + ": corr.dt must be > 0");
    if (opt_.diag_mask <= 0 || opt_.diag_mask > 7) throw std::runtime_error(std::string(tag) + ": invalid diag_mask");
    if (sel_.idx.empty()) throw std::runtime_error(std::string(tag) + ": selection is empty");
    if (opt_.remove_drift && drift_sel_.idx.empty()) throw std::runtime_error(std::string(tag) + ": drift selection is empty");
    if (opt_.corr.axis == LagAxis::TimeBin && !(opt_.corr.timebin_width > 0.0)) {
      throw std::runtime_error(std::string(tag) + ": timebin_width must be > 0 for lag_axis=timebin");
    }
    if (opt_.corr.type == "exact" && opt_.frame_end < 0) {
      throw std::runtime_error(std::string(tag) + ": correlator=exact requires finite frame_end");
    }
  }

  std::size_t exact_window_frames_() const {
    if (opt_.corr.type != "exact") return 0;
    const std::size_t w = static_cast<std::size_t>(opt_.frame_end - opt_.frame_start + 1);
    if (w < 2) throw std::runtime_error("SelfISFMeasure: exact window must have >= 2 frames");
    return w;
  }

  bool frame_in_range_(std::size_t frame_index) const {
    const std::int64_t fi = static_cast<std::int64_t>(frame_index);
    if (fi < opt_.frame_start) return false;
    if (opt_.frame_end >= 0 && fi > opt_.frame_end) return false;
    return true;
  }

  void allocate_scratch_() {
    const std::size_t nsel = sel_.idx.size();
    tmp_xx_.assign(nsel, 0.0);
    tmp_yy_.assign(nsel, 0.0);
    tmp_zz_.assign(nsel, 0.0);
    tmp_xy_.assign(nsel, 0.0);
    tmp_xz_.assign(nsel, 0.0);
    tmp_yz_.assign(nsel, 0.0);
  }

  void fill_centered_positions_(const Frame& frame) {
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
      tmp_xx_[p] = xu[i] - comx;
      tmp_yy_[p] = yu[i] - comy;
      tmp_zz_[p] = zu[i] - comz;
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

  void write_common_header_(std::ostream& ofs,
                            const char* title,
                            const char* family_name,
                            const std::string& family_values) const {
    ofs << "# PILOTS: " << title << "\n";
    ofs << "# correlator: " << opt_.corr.type << "\n";
    ofs << "# lag_axis: " << lag_axis_name(opt_.corr.axis) << "\n";
    if (opt_.corr.axis == LagAxis::TimeBin) {
      ofs << "# timebin_width: " << std::setprecision(17) << opt_.corr.timebin_width << "\n";
      ofs << "# min_pairs_per_bin: " << opt_.corr.min_pairs_per_bin << "\n";
      ofs << "# bin_merge: " << (opt_.corr.bin_merge ? 1 : 0) << "\n";
    }
    ofs << "# dt: " << std::setprecision(17) << opt_.corr.dt << " (seconds per timestep)\n";
    ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
    ofs << "# remove_drift: " << (opt_.remove_drift ? 1 : 0)
        << " (drift_group=" << drift_sel_.name << ", n=" << drift_sel_.idx.size() << ")\n";
    ofs << "# frame_range: [" << opt_.frame_start << ", " << opt_.frame_end << "]\n";
    ofs << "# components_mask: " << opt_.diag_mask << " (xx=1,yy=2,zz=4)\n";
    ofs << "# dimensions: " << count_diag_dims(opt_.diag_mask) << "\n";
    ofs << "# " << family_name << ": " << family_values << "\n";
    if (opt_.corr.type == "exact") {
      ofs << "# exact_lag_stride: " << opt_.corr.lag_stride << "\n";
    } else if (opt_.corr.type == "multitau") {
      ofs << "# multitau_channels: " << opt_.corr.mt_channels << ", levels: " << opt_.corr.mt_levels << "\n";
    }
    ofs << "# block_size: " << opt_.corr.block_size << " (0 disables SEM)\n";
  }
};

class Chi4OverlapMeasure final : public IMeasure {
public:
  Chi4OverlapMeasure(std::string instance_name,
                     std::string output_path,
                     SelectionView sel,
                     SelectionView drift_sel,
                     FamilyCommonOptions opt,
                     std::vector<double> a_values)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)),
        a_values_(std::move(a_values)) {
    sel_name_owned_ = std::string(sel.name);
    drift_name_owned_ = std::string(drift_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    drift_sel_ = SelectionView{drift_name_owned_, drift_sel.idx};

    validate_common_("Chi4OverlapMeasure");
    for (const double a : a_values_) {
      if (!(a > 0.0)) throw std::runtime_error("Chi4OverlapMeasure: each a must be > 0");
      corr_.push_back(make_correlator_t6_runtime_auto(
          sel_.idx.size(), exact_window_frames_(), opt_.corr, Chi4OverlapPairOpT6{a, opt_.diag_mask}));
    }
  }

  std::string type() const override { return "chi4_overlap"; }
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
    od.columns = {"lag", "time", "a", "q", "q2", "chi4", "count_pairs", "sem_q", "sem_q2", "n_blocks", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));

    fill_common_params_(md);
    md.params["a_values"] = join_doubles(a_values_);
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    allocate_scratch_();
    for (auto& c : corr_) c->start();
    started_ = true;
    if (!opt_.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range_(frame_index)) return;
    fill_centered_positions_(frame);
    for (auto& c : corr_) {
      c->push(frame.timestep,
              std::span<const double>(tmp_xx_.data(), tmp_xx_.size()),
              std::span<const double>(tmp_yy_.data(), tmp_yy_.size()),
              std::span<const double>(tmp_zz_.data(), tmp_zz_.size()),
              std::span<const double>(tmp_xy_.data(), tmp_xy_.size()),
              std::span<const double>(tmp_xz_.data(), tmp_xz_.size()),
              std::span<const double>(tmp_yz_.data(), tmp_yz_.size()));
    }
  }

  void flush_partial() override {
    if (opt_.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      write_common_header_(ofs, "four-point susceptibility from overlap window", "a_values", join_doubles(a_values_));
      ofs << "# columns: lag  time  a  q  q2  chi4  count_pairs  sem_q  sem_q2  n_blocks  mean_dtimestep\n";
      const double nsel = static_cast<double>(sel_.idx.size());
      for (std::size_t ia = 0; ia < a_values_.size(); ++ia) {
        const CorrelationSeriesT6 series = corr_[ia]->snapshot();
        for (std::size_t i = 0; i < series.lag.size(); ++i) {
          const double q = series.value_xx[i];
          const double q2 = series.value_yy[i];
          const double chi4 = nsel * (q2 - q * q);
          ofs << std::setprecision(17) << series.lag[i] << ' '
              << std::setprecision(17) << series.time[i] << ' '
              << std::setprecision(17) << a_values_[ia] << ' '
              << std::setprecision(17) << q << ' '
              << std::setprecision(17) << q2 << ' '
              << std::setprecision(17) << chi4 << ' '
              << series.count_pairs[i] << ' '
              << std::setprecision(17) << series.sem_xx[i] << ' '
              << std::setprecision(17) << series.sem_yy[i] << ' '
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
  FamilyCommonOptions opt_;
  std::vector<double> a_values_;
  std::vector<std::unique_ptr<ICorrelatorT6>> corr_;
  bool started_ = false;

  std::vector<double> tmp_xx_, tmp_yy_, tmp_zz_, tmp_xy_, tmp_xz_, tmp_yz_;

  void validate_common_(const char* tag) const {
    if (opt_.frame_start < 0) throw std::runtime_error(std::string(tag) + ": frame_start must be >= 0");
    if (opt_.frame_end >= 0 && opt_.frame_end < opt_.frame_start) {
      throw std::runtime_error(std::string(tag) + ": frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.corr.dt > 0.0)) throw std::runtime_error(std::string(tag) + ": corr.dt must be > 0");
    if (opt_.diag_mask <= 0 || opt_.diag_mask > 7) throw std::runtime_error(std::string(tag) + ": invalid diag_mask");
    if (sel_.idx.empty()) throw std::runtime_error(std::string(tag) + ": selection is empty");
    if (opt_.remove_drift && drift_sel_.idx.empty()) throw std::runtime_error(std::string(tag) + ": drift selection is empty");
    if (opt_.corr.axis == LagAxis::TimeBin && !(opt_.corr.timebin_width > 0.0)) {
      throw std::runtime_error(std::string(tag) + ": timebin_width must be > 0 for lag_axis=timebin");
    }
    if (opt_.corr.type == "exact" && opt_.frame_end < 0) {
      throw std::runtime_error(std::string(tag) + ": correlator=exact requires finite frame_end");
    }
  }

  std::size_t exact_window_frames_() const {
    if (opt_.corr.type != "exact") return 0;
    const std::size_t w = static_cast<std::size_t>(opt_.frame_end - opt_.frame_start + 1);
    if (w < 2) throw std::runtime_error("Chi4OverlapMeasure: exact window must have >= 2 frames");
    return w;
  }

  bool frame_in_range_(std::size_t frame_index) const {
    const std::int64_t fi = static_cast<std::int64_t>(frame_index);
    if (fi < opt_.frame_start) return false;
    if (opt_.frame_end >= 0 && fi > opt_.frame_end) return false;
    return true;
  }

  void allocate_scratch_() {
    const std::size_t nsel = sel_.idx.size();
    tmp_xx_.assign(nsel, 0.0);
    tmp_yy_.assign(nsel, 0.0);
    tmp_zz_.assign(nsel, 0.0);
    tmp_xy_.assign(nsel, 0.0);
    tmp_xz_.assign(nsel, 0.0);
    tmp_yz_.assign(nsel, 0.0);
  }

  void fill_centered_positions_(const Frame& frame) {
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
      tmp_xx_[p] = xu[i] - comx;
      tmp_yy_[p] = yu[i] - comy;
      tmp_zz_[p] = zu[i] - comz;
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

  void write_common_header_(std::ostream& ofs,
                            const char* title,
                            const char* family_name,
                            const std::string& family_values) const {
    ofs << "# PILOTS: " << title << "\n";
    ofs << "# correlator: " << opt_.corr.type << "\n";
    ofs << "# lag_axis: " << lag_axis_name(opt_.corr.axis) << "\n";
    if (opt_.corr.axis == LagAxis::TimeBin) {
      ofs << "# timebin_width: " << std::setprecision(17) << opt_.corr.timebin_width << "\n";
      ofs << "# min_pairs_per_bin: " << opt_.corr.min_pairs_per_bin << "\n";
      ofs << "# bin_merge: " << (opt_.corr.bin_merge ? 1 : 0) << "\n";
    }
    ofs << "# dt: " << std::setprecision(17) << opt_.corr.dt << " (seconds per timestep)\n";
    ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
    ofs << "# remove_drift: " << (opt_.remove_drift ? 1 : 0)
        << " (drift_group=" << drift_sel_.name << ", n=" << drift_sel_.idx.size() << ")\n";
    ofs << "# frame_range: [" << opt_.frame_start << ", " << opt_.frame_end << "]\n";
    ofs << "# components_mask: " << opt_.diag_mask << " (xx=1,yy=2,zz=4)\n";
    ofs << "# dimensions: " << count_diag_dims(opt_.diag_mask) << "\n";
    ofs << "# " << family_name << ": " << family_values << "\n";
    if (opt_.corr.type == "exact") {
      ofs << "# exact_lag_stride: " << opt_.corr.lag_stride << "\n";
    } else if (opt_.corr.type == "multitau") {
      ofs << "# multitau_channels: " << opt_.corr.mt_channels << ", levels: " << opt_.corr.mt_levels << "\n";
    }
    ofs << "# block_size: " << opt_.corr.block_size << " (0 disables SEM)\n";
  }
};

class Chi4ISFMeasure final : public IMeasure {
public:
  Chi4ISFMeasure(std::string instance_name,
                 std::string output_path,
                 SelectionView sel,
                 SelectionView drift_sel,
                 FamilyCommonOptions opt,
                 std::vector<double> q_values)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)),
        q_values_(std::move(q_values)) {
    sel_name_owned_ = std::string(sel.name);
    drift_name_owned_ = std::string(drift_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    drift_sel_ = SelectionView{drift_name_owned_, drift_sel.idx};

    validate_common_("Chi4ISFMeasure");
    for (const double q : q_values_) {
      if (!(q > 0.0)) throw std::runtime_error("Chi4ISFMeasure: each q must be > 0");
      corr_.push_back(make_correlator_t6_runtime_auto(
          sel_.idx.size(), exact_window_frames_(), opt_.corr, Chi4ISFPairOpT6{q, opt_.diag_mask}));
    }
  }

  std::string type() const override { return "chi4_isf"; }
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
    od.columns = {"lag", "time", "q", "fs", "fs2", "chi4", "count_pairs", "sem_fs", "sem_fs2", "n_blocks", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));

    fill_common_params_(md);
    md.params["q_values"] = join_doubles(q_values_);
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    allocate_scratch_();
    for (auto& c : corr_) c->start();
    started_ = true;
    if (!opt_.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range_(frame_index)) return;
    fill_centered_positions_(frame);
    for (auto& c : corr_) {
      c->push(frame.timestep,
              std::span<const double>(tmp_xx_.data(), tmp_xx_.size()),
              std::span<const double>(tmp_yy_.data(), tmp_yy_.size()),
              std::span<const double>(tmp_zz_.data(), tmp_zz_.size()),
              std::span<const double>(tmp_xy_.data(), tmp_xy_.size()),
              std::span<const double>(tmp_xz_.data(), tmp_xz_.size()),
              std::span<const double>(tmp_yz_.data(), tmp_yz_.size()));
    }
  }

  void flush_partial() override {
    if (opt_.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      write_common_header_(ofs, "four-point susceptibility from self-ISF fluctuations", "q_values", join_doubles(q_values_));
      ofs << "# columns: lag  time  q  fs  fs2  chi4  count_pairs  sem_fs  sem_fs2  n_blocks  mean_dtimestep\n";
      const double nsel = static_cast<double>(sel_.idx.size());
      for (std::size_t iq = 0; iq < q_values_.size(); ++iq) {
        const CorrelationSeriesT6 series = corr_[iq]->snapshot();
        for (std::size_t i = 0; i < series.lag.size(); ++i) {
          const double fs = series.value_xx[i];
          const double fs2 = series.value_yy[i];
          const double chi4 = nsel * (fs2 - fs * fs);
          ofs << std::setprecision(17) << series.lag[i] << ' '
              << std::setprecision(17) << series.time[i] << ' '
              << std::setprecision(17) << q_values_[iq] << ' '
              << std::setprecision(17) << fs << ' '
              << std::setprecision(17) << fs2 << ' '
              << std::setprecision(17) << chi4 << ' '
              << series.count_pairs[i] << ' '
              << std::setprecision(17) << series.sem_xx[i] << ' '
              << std::setprecision(17) << series.sem_yy[i] << ' '
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
  FamilyCommonOptions opt_;
  std::vector<double> q_values_;
  std::vector<std::unique_ptr<ICorrelatorT6>> corr_;
  bool started_ = false;

  std::vector<double> tmp_xx_, tmp_yy_, tmp_zz_, tmp_xy_, tmp_xz_, tmp_yz_;

  void validate_common_(const char* tag) const {
    if (opt_.frame_start < 0) throw std::runtime_error(std::string(tag) + ": frame_start must be >= 0");
    if (opt_.frame_end >= 0 && opt_.frame_end < opt_.frame_start) {
      throw std::runtime_error(std::string(tag) + ": frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.corr.dt > 0.0)) throw std::runtime_error(std::string(tag) + ": corr.dt must be > 0");
    if (opt_.diag_mask <= 0 || opt_.diag_mask > 7) throw std::runtime_error(std::string(tag) + ": invalid diag_mask");
    if (sel_.idx.empty()) throw std::runtime_error(std::string(tag) + ": selection is empty");
    if (opt_.remove_drift && drift_sel_.idx.empty()) throw std::runtime_error(std::string(tag) + ": drift selection is empty");
    if (opt_.corr.axis == LagAxis::TimeBin && !(opt_.corr.timebin_width > 0.0)) {
      throw std::runtime_error(std::string(tag) + ": timebin_width must be > 0 for lag_axis=timebin");
    }
    if (opt_.corr.type == "exact" && opt_.frame_end < 0) {
      throw std::runtime_error(std::string(tag) + ": correlator=exact requires finite frame_end");
    }
  }

  std::size_t exact_window_frames_() const {
    if (opt_.corr.type != "exact") return 0;
    const std::size_t w = static_cast<std::size_t>(opt_.frame_end - opt_.frame_start + 1);
    if (w < 2) throw std::runtime_error("Chi4ISFMeasure: exact window must have >= 2 frames");
    return w;
  }

  bool frame_in_range_(std::size_t frame_index) const {
    const std::int64_t fi = static_cast<std::int64_t>(frame_index);
    if (fi < opt_.frame_start) return false;
    if (opt_.frame_end >= 0 && fi > opt_.frame_end) return false;
    return true;
  }

  void allocate_scratch_() {
    const std::size_t nsel = sel_.idx.size();
    tmp_xx_.assign(nsel, 0.0);
    tmp_yy_.assign(nsel, 0.0);
    tmp_zz_.assign(nsel, 0.0);
    tmp_xy_.assign(nsel, 0.0);
    tmp_xz_.assign(nsel, 0.0);
    tmp_yz_.assign(nsel, 0.0);
  }

  void fill_centered_positions_(const Frame& frame) {
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
      tmp_xx_[p] = xu[i] - comx;
      tmp_yy_[p] = yu[i] - comy;
      tmp_zz_[p] = zu[i] - comz;
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

  void write_common_header_(std::ostream& ofs,
                            const char* title,
                            const char* family_name,
                            const std::string& family_values) const {
    ofs << "# PILOTS: " << title << "\n";
    ofs << "# correlator: " << opt_.corr.type << "\n";
    ofs << "# lag_axis: " << lag_axis_name(opt_.corr.axis) << "\n";
    if (opt_.corr.axis == LagAxis::TimeBin) {
      ofs << "# timebin_width: " << std::setprecision(17) << opt_.corr.timebin_width << "\n";
      ofs << "# min_pairs_per_bin: " << opt_.corr.min_pairs_per_bin << "\n";
      ofs << "# bin_merge: " << (opt_.corr.bin_merge ? 1 : 0) << "\n";
    }
    ofs << "# dt: " << std::setprecision(17) << opt_.corr.dt << " (seconds per timestep)\n";
    ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
    ofs << "# remove_drift: " << (opt_.remove_drift ? 1 : 0)
        << " (drift_group=" << drift_sel_.name << ", n=" << drift_sel_.idx.size() << ")\n";
    ofs << "# frame_range: [" << opt_.frame_start << ", " << opt_.frame_end << "]\n";
    ofs << "# components_mask: " << opt_.diag_mask << " (xx=1,yy=2,zz=4)\n";
    ofs << "# dimensions: " << count_diag_dims(opt_.diag_mask) << "\n";
    ofs << "# " << family_name << ": " << family_values << "\n";
    if (opt_.corr.type == "exact") {
      ofs << "# exact_lag_stride: " << opt_.corr.lag_stride << "\n";
    } else if (opt_.corr.type == "multitau") {
      ofs << "# multitau_channels: " << opt_.corr.mt_channels << ", levels: " << opt_.corr.mt_levels << "\n";
    }
    ofs << "# block_size: " << opt_.corr.block_size << " (0 disables SEM)\n";
  }
};

MeasureCapabilities scattering_family_caps(const IniConfig& cfg,
                                           const std::string& section,
                                           const std::string& instance,
                                           const MeasureBuildEnv& env) {
  (void)instance;
  (void)env;
  MeasureCapabilities caps;
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

  return caps;
}

FamilyCommonOptions make_family_common_options(const IniConfig& cfg,
                                                const std::string& section,
                                                const MeasureBuildEnv& env,
                                                int& diag_mask_out) {
  const CorrelatorSpec corr_spec = parse_correlator_spec(cfg, section, env.dt);
  const std::int64_t frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  const std::int64_t frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  if (frame_start < 0) throw std::runtime_error("frame_start must be >= 0");
  if (frame_end >= 0 && frame_end < frame_start) {
    throw std::runtime_error("frame_end must be -1 or >= frame_start");
  }

  const std::string comp_s = cfg.get_string(section, "components", std::optional<std::string>("xxyyzz"));
  diag_mask_out = parse_diag_mask(comp_s);

  FamilyCommonOptions opt;
  opt.frame_start = frame_start;
  opt.frame_end = (corr_spec.type == "exact")
      ? resolve_exact_frame_end(corr_spec, env.follow, frame_start, frame_end, env.input_path)
      : frame_end;
  opt.diag_mask = diag_mask_out;
  opt.remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(true));
  opt.dry_run = env.dry_run;
  opt.corr = corr_spec;
  return opt;
}

std::pair<SelectionView, SelectionView> resolve_scattering_selections(const IniConfig& cfg,
                                                                      const std::string& section,
                                                                      const MeasureBuildEnv& env) {
  if (!env.first_frame) throw std::runtime_error("scattering family factory: first_frame is null");
  if (!env.selection_provider) throw std::runtime_error("scattering family factory: SelectionProvider is missing");
  const Frame& frame0 = *env.first_frame;
  const std::string group_ref = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo_group_ref = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string combine_expr = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(true));
  const std::string drift_group_ref = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));

  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group_ref, topo_group_ref, combine_expr, "scattering family");
  SelectionView drift_sel = remove_drift
      ? get_static_group_view(*env.selection_provider, frame0, drift_group_ref, "scattering family drift_group")
      : get_static_group_view(*env.selection_provider, frame0, "all", "scattering family drift_group");
  return {sel, drift_sel};
}

fs::path resolve_output_path(const IniConfig& cfg,
                             const std::string& section,
                             const MeasureBuildEnv& env,
                             const std::string& default_name) {
  const std::string out_file = cfg.get_string(section, "output", std::optional<std::string>(default_name));
  const std::string outdir_s = cfg.get_string(section, "output_dir", std::optional<std::string>(""));
  const fs::path output_dir = outdir_s.empty() ? env.output_dir_general : resolve_path(env.cfg_dir, outdir_s);
  if (!env.dry_run) fs::create_directories(output_dir);
  return (output_dir / out_file).lexically_normal();
}

MeasureCapabilities self_isf_caps(const IniConfig& cfg,
                                  const std::string& section,
                                  const std::string& instance,
                                  const MeasureBuildEnv& env) {
  (void)instance;
  return scattering_family_caps(cfg, section, instance, env);
}

std::unique_ptr<IMeasure> self_isf_create(const IniConfig& cfg,
                                          const std::string& section,
                                          const std::string& instance,
                                          const MeasureBuildEnv& env,
                                          const SystemContext& sysctx) {
  (void)sysctx;
  auto [sel, drift_sel] = resolve_scattering_selections(cfg, section, env);
  int diag_mask = 7;
  FamilyCommonOptions opt = make_family_common_options(cfg, section, env, diag_mask);
  const auto q_values = parse_double_list(cfg, section, "q_values");
  const fs::path out_path = resolve_output_path(cfg, section, env, "self_isf.dat");
  return std::make_unique<SelfISFMeasure>(instance, out_path.string(), sel, drift_sel, opt, q_values);
}

MeasureCapabilities chi4_overlap_caps(const IniConfig& cfg,
                                      const std::string& section,
                                      const std::string& instance,
                                      const MeasureBuildEnv& env) {
  (void)instance;
  return scattering_family_caps(cfg, section, instance, env);
}

std::unique_ptr<IMeasure> chi4_overlap_create(const IniConfig& cfg,
                                              const std::string& section,
                                              const std::string& instance,
                                              const MeasureBuildEnv& env,
                                              const SystemContext& sysctx) {
  (void)sysctx;
  auto [sel, drift_sel] = resolve_scattering_selections(cfg, section, env);
  int diag_mask = 7;
  FamilyCommonOptions opt = make_family_common_options(cfg, section, env, diag_mask);
  const auto a_values = parse_double_list(cfg, section, "a_values");
  const fs::path out_path = resolve_output_path(cfg, section, env, "chi4_overlap.dat");
  return std::make_unique<Chi4OverlapMeasure>(instance, out_path.string(), sel, drift_sel, opt, a_values);
}

MeasureCapabilities chi4_isf_caps(const IniConfig& cfg,
                                  const std::string& section,
                                  const std::string& instance,
                                  const MeasureBuildEnv& env) {
  (void)instance;
  return scattering_family_caps(cfg, section, instance, env);
}

std::unique_ptr<IMeasure> chi4_isf_create(const IniConfig& cfg,
                                          const std::string& section,
                                          const std::string& instance,
                                          const MeasureBuildEnv& env,
                                          const SystemContext& sysctx) {
  (void)sysctx;
  auto [sel, drift_sel] = resolve_scattering_selections(cfg, section, env);
  int diag_mask = 7;
  FamilyCommonOptions opt = make_family_common_options(cfg, section, env, diag_mask);
  const auto q_values = parse_double_list(cfg, section, "q_values");
  const fs::path out_path = resolve_output_path(cfg, section, env, "chi4_isf.dat");
  return std::make_unique<Chi4ISFMeasure>(instance, out_path.string(), sel, drift_sel, opt, q_values);
}

static MeasureRegistrar g_register_self_isf("self_isf", &self_isf_caps, &self_isf_create);
static MeasureRegistrar g_register_chi4_overlap("chi4_overlap", &chi4_overlap_caps, &chi4_overlap_create);
static MeasureRegistrar g_register_chi4_isf("chi4_isf", &chi4_isf_caps, &chi4_isf_create);

} // namespace
} // namespace pilots
