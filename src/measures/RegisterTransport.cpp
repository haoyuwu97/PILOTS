#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
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
using measure_ext::dfield_to_i64;
using measure_ext::dstr;
using measure_ext::get_static_combined_view;
using measure_ext::parse_diag_mask;
using measure_ext::resolve_exact_frame_end;
using measure_ext::resolve_path;
using measure_ext::x_unit_for_axis;

struct TransportRangeOptions {
  std::int64_t frame_start = 0;
  std::int64_t frame_end = -1;
  bool dry_run = false;
};

struct CurrentACFOptions {
  TransportRangeOptions range;
  CorrelatorSpec corr;
  int diag_mask = 7;
  std::string q_field = "q";
  std::string vector_x_field = "vx";
  std::string vector_y_field = "vy";
  std::string vector_z_field = "vz";
  double kbt = 1.0;
};

struct ChargeEinsteinOptions {
  TransportRangeOptions range;
  CorrelatorSpec corr;
  int diag_mask = 7;
  std::string q_field = "q";
  double kbt = 1.0;
};

struct SiteHopOptions {
  TransportRangeOptions range;
  std::string site_field = "site";
  std::int64_t unassigned_site = -1;
  std::size_t min_persistence_frames = 1;
  bool ignore_unassigned = true;
  double dt = 0.0;
};

inline void append_static_selection_caps(const IniConfig& cfg,
                                         const std::string& section,
                                         MeasureCapabilities& caps) {
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_identity_consistent = true;
  caps.scale = ScaleCompatibility{true, true, true};
  caps.group_refs.push_back(cfg.get_string(section, "group", std::optional<std::string>("all")));
}

inline void append_integer_like_field_cap(MeasureCapabilities& caps,
                                          const std::string& field) {
  if (field == "id" || field == "mol") {
    caps.requires_i64fields.push_back(field);
  } else if (field == "type") {
    caps.requires_intfields.push_back(field);
  } else {
    caps.requires_dfields.push_back(field);
  }
}

inline std::vector<std::int64_t> integer_like_field_to_i64(const Frame& frame,
                                                           const std::string& field) {
  if (field == "id" || field == "mol") {
    const auto v = frame.require_i64field(field);
    return std::vector<std::int64_t>(v.begin(), v.end());
  }
  if (field == "type") {
    const auto v = frame.require_intfield(field);
    std::vector<std::int64_t> out(v.size(), 0);
    for (std::size_t i = 0; i < v.size(); ++i) out[i] = static_cast<std::int64_t>(v[i]);
    return out;
  }
  return dfield_to_i64(frame.require_dfield(field), field, true);
}

inline std::vector<double> selected_field_copy(SelectionView sel,
                                               std::span<const double> field) {
  std::vector<double> out(sel.idx.size(), 0.0);
#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
  for (std::size_t p = 0; p < sel.idx.size(); ++p) {
    out[p] = field[sel.idx[p]];
  }
  return out;
}

inline std::vector<double> running_trapezoid(const std::vector<double>& x,
                                             const std::vector<double>& y) {
  const std::size_t n = std::min(x.size(), y.size());
  std::vector<double> out(n, 0.0);
  for (std::size_t i = 1; i < n; ++i) {
    const double dx = x[i] - x[i - 1];
    out[i] = out[i - 1] + 0.5 * dx * (y[i] + y[i - 1]);
  }
  return out;
}

inline const std::vector<double>& series_time_axis(const CorrelationSeriesT6& series) {
  if (!series.time.empty()) return series.time;
  return series.lag;
}

struct CurrentPairOpT6 {
  int diag_mask = 7;

  Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const {
    if (cur.xx.size() != 1 || cur.yy.size() != 1 || cur.zz.size() != 1 ||
        org.xx.size() != 1 || org.yy.size() != 1 || org.zz.size() != 1) {
      throw std::runtime_error("CurrentPairOpT6: expected a single current vector signal");
    }
    Tensor6 out;
    if ((diag_mask & 1) != 0) out.v[T6_XX] = cur.xx[0] * org.xx[0];
    if ((diag_mask & 2) != 0) out.v[T6_YY] = cur.yy[0] * org.yy[0];
    if ((diag_mask & 4) != 0) out.v[T6_ZZ] = cur.zz[0] * org.zz[0];
    return out;
  }
};

struct ChargeEinsteinPairOpT6 {
  std::vector<double> q;
  int diag_mask = 7;

  Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const {
    const std::size_t n = q.size();
    if (cur.xx.size() != n || cur.yy.size() != n || cur.zz.size() != n ||
        org.xx.size() != n || org.yy.size() != n || org.zz.size() != n) {
      throw std::runtime_error("ChargeEinsteinPairOpT6: slot sizes mismatch with charge vector");
    }

    double mx = 0.0, my = 0.0, mz = 0.0;
    double sx = 0.0, sy = 0.0, sz = 0.0;
#if PILOTS_HAS_OPENMP
#pragma omp parallel for reduction(+:mx,my,mz,sx,sy,sz)
#endif
    for (std::size_t p = 0; p < n; ++p) {
      const double dx = cur.xx[p] - org.xx[p];
      const double dy = cur.yy[p] - org.yy[p];
      const double dz = cur.zz[p] - org.zz[p];
      const double qi = q[p];
      const double qi2 = qi * qi;
      if ((diag_mask & 1) != 0) {
        mx += qi * dx;
        sx += qi2 * dx * dx;
      }
      if ((diag_mask & 2) != 0) {
        my += qi * dy;
        sy += qi2 * dy * dy;
      }
      if ((diag_mask & 4) != 0) {
        mz += qi * dz;
        sz += qi2 * dz * dz;
      }
    }

    Tensor6 out;
    out.v[T6_XX] = mx * mx;
    out.v[T6_YY] = my * my;
    out.v[T6_ZZ] = mz * mz;
    out.v[T6_XY] = sx;
    out.v[T6_XZ] = sy;
    out.v[T6_YZ] = sz;
    return out;
  }
};

class CurrentACFMeasure final : public IMeasure {
public:
  CurrentACFMeasure(std::string instance_name,
                    std::string output_path,
                    SelectionView sel,
                    CurrentACFOptions opt)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    validate_();
    corr_ = make_correlator_t6_runtime_auto(1, exact_window_frames_(), opt_.corr, CurrentPairOpT6{opt_.diag_mask});
  }

  std::string type() const override { return "current_acf"; }
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
    od.x_axis = pilots::lag_axis_name(opt_.corr.axis);
    od.x_unit = x_unit_for_axis(opt_.corr.axis);
    od.columns = {
        "lag", "time",
        "cjj_xx", "cjj_yy", "cjj_zz", "cjj_trace",
        "sigma_xx", "sigma_yy", "sigma_zz", "sigma_gk",
        "count_pairs", "n_blocks", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));

    md.params["components"] = components_name_();
    md.params["q_field"] = opt_.q_field;
    md.params["vector_x_field"] = opt_.vector_x_field;
    md.params["vector_y_field"] = opt_.vector_y_field;
    md.params["vector_z_field"] = opt_.vector_z_field;
    md.params["kbt"] = dstr(opt_.kbt);
    md.params["frame_start"] = std::to_string(opt_.range.frame_start);
    md.params["frame_end"] = std::to_string(opt_.range.frame_end);
    md.params["correlator"] = opt_.corr.type;
    md.params["lag_axis"] = pilots::lag_axis_name(opt_.corr.axis);
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    jx_.assign(1, 0.0);
    jy_.assign(1, 0.0);
    jz_.assign(1, 0.0);
    zero_.assign(1, 0.0);
    corr_->start();
    started_ = true;
    if (!opt_.range.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range_(frame_index)) return;

    const auto q = frame.require_dfield(opt_.q_field);
    const auto ux = frame.require_dfield(opt_.vector_x_field);
    const auto uy = frame.require_dfield(opt_.vector_y_field);
    const auto uz = frame.require_dfield(opt_.vector_z_field);

    double sumx = 0.0, sumy = 0.0, sumz = 0.0;
#if PILOTS_HAS_OPENMP
#pragma omp parallel for reduction(+:sumx,sumy,sumz)
#endif
    for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
      const std::size_t i = sel_.idx[p];
      const double qi = q[i];
      sumx += qi * ux[i];
      sumy += qi * uy[i];
      sumz += qi * uz[i];
    }
    jx_[0] = sumx;
    jy_[0] = sumy;
    jz_[0] = sumz;
    corr_->push(frame.timestep,
                std::span<const double>(jx_.data(), 1),
                std::span<const double>(jy_.data(), 1),
                std::span<const double>(jz_.data(), 1),
                std::span<const double>(zero_.data(), 1),
                std::span<const double>(zero_.data(), 1),
                std::span<const double>(zero_.data(), 1));
    volume_sum_ += frame.box.lx() * frame.box.ly() * frame.box.lz();
    ++n_frames_used_;
  }

  void flush_partial() override {
    write_output_();
  }

  void finalize() override {
    write_output_();
  }

private:
  void validate_() const {
    if (sel_.idx.empty()) throw std::runtime_error("current_acf: selection is empty");
    if (opt_.range.frame_start < 0) throw std::runtime_error("current_acf: frame_start must be >= 0");
    if (opt_.range.frame_end >= 0 && opt_.range.frame_end < opt_.range.frame_start) {
      throw std::runtime_error("current_acf: frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.kbt > 0.0)) throw std::runtime_error("current_acf: kbt must be > 0");
    if (count_diag_dims(opt_.diag_mask) <= 0) {
      throw std::runtime_error("current_acf: at least one component must be selected");
    }
  }

  std::size_t exact_window_frames_() const {
    if (opt_.corr.type != "exact") return 0;
    return static_cast<std::size_t>(opt_.range.frame_end - opt_.range.frame_start + 1);
  }

  bool frame_in_range_(std::size_t frame_index) const {
    const std::int64_t fi = static_cast<std::int64_t>(frame_index);
    if (fi < opt_.range.frame_start) return false;
    if (opt_.range.frame_end >= 0 && fi > opt_.range.frame_end) return false;
    return true;
  }

  std::string components_name_() const {
    switch (opt_.diag_mask) {
      case 1: return "xx";
      case 2: return "yy";
      case 4: return "zz";
      case 3: return "xxyy";
      case 5: return "xxzz";
      case 6: return "yyzz";
      default: return "xxyyzz";
    }
  }

  void write_output_() const {
    const CorrelationSeriesT6 series = corr_->snapshot();
    const auto& t = series_time_axis(series);
    const std::vector<double> int_x = running_trapezoid(t, series.value_xx);
    const std::vector<double> int_y = running_trapezoid(t, series.value_yy);
    const std::vector<double> int_z = running_trapezoid(t, series.value_zz);
    const double mean_volume = (n_frames_used_ > 0) ? (volume_sum_ / static_cast<double>(n_frames_used_)) : 0.0;
    const double pref = (mean_volume > 0.0) ? (1.0 / (mean_volume * opt_.kbt)) : 0.0;
    const double dinv = 1.0 / static_cast<double>(count_diag_dims(opt_.diag_mask));

    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << std::setprecision(17);
      write_header_(ofs);
      const std::size_t n = series.lag.size();
      for (std::size_t i = 0; i < n; ++i) {
        const double cxx = (i < series.value_xx.size()) ? series.value_xx[i] : 0.0;
        const double cyy = (i < series.value_yy.size()) ? series.value_yy[i] : 0.0;
        const double czz = (i < series.value_zz.size()) ? series.value_zz[i] : 0.0;
        const double sx = (i < int_x.size()) ? pref * int_x[i] : 0.0;
        const double sy = (i < int_y.size()) ? pref * int_y[i] : 0.0;
        const double sz = (i < int_z.size()) ? pref * int_z[i] : 0.0;
        const double ctrace = cxx + cyy + czz;
        const double sigma = dinv * (sx + sy + sz);
        ofs << series.lag[i] << ' '
            << ((i < t.size()) ? t[i] : 0.0) << ' '
            << cxx << ' ' << cyy << ' ' << czz << ' ' << ctrace << ' '
            << sx << ' ' << sy << ' ' << sz << ' ' << sigma << ' '
            << ((i < series.count_pairs.size()) ? series.count_pairs[i] : 0) << ' '
            << ((i < series.n_blocks.size()) ? series.n_blocks[i] : 0) << ' '
            << ((i < series.mean_dtimestep.size()) ? series.mean_dtimestep[i] : 0.0) << '\n';
      }
    });
  }

  void write_header_(std::ostream& ofs) const {
    ofs << "# PILOTS: current autocorrelation / Green-Kubo conductivity\n";
    ofs << "# selection: " << sel_.name << "\n";
    ofs << "# components: " << components_name_() << "\n";
    ofs << "# q_field: " << opt_.q_field
        << ", vector_fields: (" << opt_.vector_x_field << ',' << opt_.vector_y_field << ',' << opt_.vector_z_field << ")\n";
    ofs << "# correlator: " << opt_.corr.type << ", lag_axis: " << pilots::lag_axis_name(opt_.corr.axis) << "\n";
    ofs << "# kbt: " << opt_.kbt << "\n";
    ofs << "# columns: lag  time  cjj_xx  cjj_yy  cjj_zz  cjj_trace  sigma_xx  sigma_yy  sigma_zz  sigma_gk  count_pairs  n_blocks  mean_dtimestep\n";
  }

  std::string instance_name_;
  std::string output_path_;
  CurrentACFOptions opt_;
  std::string sel_name_owned_;
  SelectionView sel_{};
  mutable std::unique_ptr<ICorrelatorT6> corr_;
  std::vector<double> jx_, jy_, jz_, zero_;
  bool started_ = false;
  double volume_sum_ = 0.0;
  std::size_t n_frames_used_ = 0;
};

class ChargeEinsteinMeasure final : public IMeasure {
public:
  ChargeEinsteinMeasure(std::string instance_name,
                        std::string output_path,
                        SelectionView sel,
                        ChargeEinsteinOptions opt,
                        std::vector<double> q_selected)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)),
        q_(std::move(q_selected)) {
    sel_name_owned_ = std::string(sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    validate_();
    corr_ = make_correlator_t6_runtime_auto(sel_.idx.size(), exact_window_frames_(), opt_.corr,
                                            ChargeEinsteinPairOpT6{q_, opt_.diag_mask});
  }

  std::string type() const override { return "conductivity_einstein"; }
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
    od.x_axis = pilots::lag_axis_name(opt_.corr.axis);
    od.x_unit = x_unit_for_axis(opt_.corr.axis);
    od.columns = {
        "lag", "time",
        "m_xx", "m_yy", "m_zz", "m_total",
        "self_xx", "self_yy", "self_zz", "self_total", "distinct_total",
        "sigma_xx", "sigma_yy", "sigma_zz", "sigma_total", "sigma_self", "sigma_distinct", "haven_ratio",
        "count_pairs", "sem_m_xx", "sem_m_yy", "sem_m_zz", "n_blocks", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));

    md.params["components"] = components_name_();
    md.params["q_field"] = opt_.q_field;
    md.params["kbt"] = dstr(opt_.kbt);
    md.params["frame_start"] = std::to_string(opt_.range.frame_start);
    md.params["frame_end"] = std::to_string(opt_.range.frame_end);
    md.params["correlator"] = opt_.corr.type;
    md.params["lag_axis"] = pilots::lag_axis_name(opt_.corr.axis);
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    x_.assign(sel_.idx.size(), 0.0);
    y_.assign(sel_.idx.size(), 0.0);
    z_.assign(sel_.idx.size(), 0.0);
    zero_.assign(sel_.idx.size(), 0.0);
    corr_->start();
    if (!opt_.range.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range_(frame_index)) return;
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");
#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
    for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
      const std::size_t i = sel_.idx[p];
      x_[p] = xu[i];
      y_[p] = yu[i];
      z_[p] = zu[i];
    }
    corr_->push(frame.timestep,
                std::span<const double>(x_.data(), x_.size()),
                std::span<const double>(y_.data(), y_.size()),
                std::span<const double>(z_.data(), z_.size()),
                std::span<const double>(zero_.data(), zero_.size()),
                std::span<const double>(zero_.data(), zero_.size()),
                std::span<const double>(zero_.data(), zero_.size()));
    volume_sum_ += frame.box.lx() * frame.box.ly() * frame.box.lz();
    ++n_frames_used_;
  }

  void flush_partial() override {
    write_output_();
  }

  void finalize() override {
    write_output_();
  }

private:
  void validate_() const {
    if (sel_.idx.empty()) throw std::runtime_error("conductivity_einstein: selection is empty");
    if (q_.size() != sel_.idx.size()) {
      throw std::runtime_error("conductivity_einstein: selected charge vector size does not match selection size");
    }
    if (opt_.range.frame_start < 0) throw std::runtime_error("conductivity_einstein: frame_start must be >= 0");
    if (opt_.range.frame_end >= 0 && opt_.range.frame_end < opt_.range.frame_start) {
      throw std::runtime_error("conductivity_einstein: frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.kbt > 0.0)) throw std::runtime_error("conductivity_einstein: kbt must be > 0");
    if (count_diag_dims(opt_.diag_mask) <= 0) {
      throw std::runtime_error("conductivity_einstein: at least one component must be selected");
    }
  }

  std::size_t exact_window_frames_() const {
    if (opt_.corr.type != "exact") return 0;
    return static_cast<std::size_t>(opt_.range.frame_end - opt_.range.frame_start + 1);
  }

  bool frame_in_range_(std::size_t frame_index) const {
    const std::int64_t fi = static_cast<std::int64_t>(frame_index);
    if (fi < opt_.range.frame_start) return false;
    if (opt_.range.frame_end >= 0 && fi > opt_.range.frame_end) return false;
    return true;
  }

  std::string components_name_() const {
    switch (opt_.diag_mask) {
      case 1: return "xx";
      case 2: return "yy";
      case 4: return "zz";
      case 3: return "xxyy";
      case 5: return "xxzz";
      case 6: return "yyzz";
      default: return "xxyyzz";
    }
  }

  void write_output_() const {
    const CorrelationSeriesT6 series = corr_->snapshot();
    const auto& t = series_time_axis(series);
    const double mean_volume = (n_frames_used_ > 0) ? (volume_sum_ / static_cast<double>(n_frames_used_)) : 0.0;
    const double pref = (mean_volume > 0.0) ? (1.0 / (mean_volume * opt_.kbt)) : 0.0;
    const double dinv = 1.0 / static_cast<double>(count_diag_dims(opt_.diag_mask));

    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << std::setprecision(17);
      write_header_(ofs);
      const std::size_t n = series.lag.size();
      for (std::size_t i = 0; i < n; ++i) {
        const double tau = (i < t.size()) ? t[i] : 0.0;
        const double mxx = (i < series.value_xx.size()) ? series.value_xx[i] : 0.0;
        const double myy = (i < series.value_yy.size()) ? series.value_yy[i] : 0.0;
        const double mzz = (i < series.value_zz.size()) ? series.value_zz[i] : 0.0;
        const double sxx = (i < series.value_xy.size()) ? series.value_xy[i] : 0.0;
        const double syy = (i < series.value_xz.size()) ? series.value_xz[i] : 0.0;
        const double szz = (i < series.value_yz.size()) ? series.value_yz[i] : 0.0;
        const double mtotal = mxx + myy + mzz;
        const double stotal = sxx + syy + szz;
        const double distinct = mtotal - stotal;
        double sigma_xx = 0.0, sigma_yy = 0.0, sigma_zz = 0.0;
        double sigma_total = 0.0, sigma_self = 0.0, sigma_distinct = 0.0;
        if (tau > 0.0 && pref > 0.0) {
          sigma_xx = pref * mxx / (2.0 * tau);
          sigma_yy = pref * myy / (2.0 * tau);
          sigma_zz = pref * mzz / (2.0 * tau);
          sigma_total = dinv * (sigma_xx + sigma_yy + sigma_zz);
          sigma_self = pref * stotal / (2.0 * static_cast<double>(count_diag_dims(opt_.diag_mask)) * tau);
          sigma_distinct = sigma_total - sigma_self;
        }
        const double haven = (sigma_self > 0.0) ? (sigma_total / sigma_self) : 0.0;
        ofs << series.lag[i] << ' '
            << tau << ' '
            << mxx << ' ' << myy << ' ' << mzz << ' ' << mtotal << ' '
            << sxx << ' ' << syy << ' ' << szz << ' ' << stotal << ' ' << distinct << ' '
            << sigma_xx << ' ' << sigma_yy << ' ' << sigma_zz << ' '
            << sigma_total << ' ' << sigma_self << ' ' << sigma_distinct << ' ' << haven << ' '
            << ((i < series.count_pairs.size()) ? series.count_pairs[i] : 0) << ' '
            << ((i < series.sem_xx.size()) ? series.sem_xx[i] : 0.0) << ' '
            << ((i < series.sem_yy.size()) ? series.sem_yy[i] : 0.0) << ' '
            << ((i < series.sem_zz.size()) ? series.sem_zz[i] : 0.0) << ' '
            << ((i < series.n_blocks.size()) ? series.n_blocks[i] : 0) << ' '
            << ((i < series.mean_dtimestep.size()) ? series.mean_dtimestep[i] : 0.0) << '\n';
      }
    });
  }

  void write_header_(std::ostream& ofs) const {
    ofs << "# PILOTS: collective-charge Einstein conductivity / Nernst-Einstein decomposition\n";
    ofs << "# selection: " << sel_.name << "\n";
    ofs << "# components: " << components_name_() << "\n";
    ofs << "# q_field: " << opt_.q_field << "\n";
    ofs << "# correlator: " << opt_.corr.type << ", lag_axis: " << pilots::lag_axis_name(opt_.corr.axis) << "\n";
    ofs << "# kbt: " << opt_.kbt << "\n";
    ofs << "# columns: lag  time  m_xx  m_yy  m_zz  m_total  self_xx  self_yy  self_zz  self_total  distinct_total  sigma_xx  sigma_yy  sigma_zz  sigma_total  sigma_self  sigma_distinct  haven_ratio  count_pairs  sem_m_xx  sem_m_yy  sem_m_zz  n_blocks  mean_dtimestep\n";
  }

  std::string instance_name_;
  std::string output_path_;
  ChargeEinsteinOptions opt_;
  std::vector<double> q_;
  std::string sel_name_owned_;
  SelectionView sel_{};
  mutable std::unique_ptr<ICorrelatorT6> corr_;
  std::vector<double> x_, y_, z_, zero_;
  double volume_sum_ = 0.0;
  std::size_t n_frames_used_ = 0;
};

class SiteHopStatsMeasure final : public IMeasure {
public:
  SiteHopStatsMeasure(std::string instance_name,
                      std::string summary_path,
                      std::string transitions_path,
                      std::string dwell_path,
                      SelectionView sel,
                      SiteHopOptions opt)
      : instance_name_(std::move(instance_name)),
        summary_path_(std::move(summary_path)),
        transitions_path_(std::move(transitions_path)),
        dwell_path_(std::move(dwell_path)),
        opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    validate_();
  }

  std::string type() const override { return "site_hop_stats"; }
  std::string instance_name() const override { return instance_name_; }

  output::MeasureDescriptor describe() const override {
    output::MeasureDescriptor md;
    md.instance = instance_name_;
    md.type = type();
    md.selection = std::string(sel_.name);
    md.n_selected = sel_.idx.size();

    output::OutputFileDescriptor od0;
    od0.path = summary_path_;
    od0.format = "text";
    od0.x_axis = "none";
    od0.x_unit = "none";
    od0.columns = {
        "n_particles", "n_frames", "n_hops_total", "hop_rate_per_particle_per_frame",
        "hop_rate_per_particle_per_time", "mean_dwell_frames", "mean_dwell_time"};
    md.outputs.push_back(std::move(od0));

    output::OutputFileDescriptor od1;
    od1.path = transitions_path_;
    od1.format = "text";
    od1.x_axis = "site_from";
    od1.x_unit = "id";
    od1.columns = {"site_from", "site_to", "count", "probability"};
    md.outputs.push_back(std::move(od1));

    output::OutputFileDescriptor od2;
    od2.path = dwell_path_;
    od2.format = "text";
    od2.x_axis = "dwell_frames";
    od2.x_unit = "frames";
    od2.columns = {"dwell_frames", "dwell_time", "count", "probability"};
    md.outputs.push_back(std::move(od2));

    md.params["site_field"] = opt_.site_field;
    md.params["unassigned_site"] = std::to_string(opt_.unassigned_site);
    md.params["min_persistence_frames"] = std::to_string(opt_.min_persistence_frames);
    md.params["ignore_unassigned"] = opt_.ignore_unassigned ? "1" : "0";
    md.params["frame_start"] = std::to_string(opt_.range.frame_start);
    md.params["frame_end"] = std::to_string(opt_.range.frame_end);
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    states_.assign(sel_.idx.size(), ParticleState{});
    if (!opt_.range.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range_(frame_index)) return;
    const auto site_all = integer_like_field_to_i64(frame, opt_.site_field);
    for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
      const std::size_t i = sel_.idx[p];
      const std::int64_t obs = site_all[i];
      auto& st = states_[p];

      if (!st.initialized) {
        st.initialized = true;
        st.committed_site = obs;
        st.entry_frame = frame_index;
        st.entry_timestep = frame.timestep;
        st.candidate_site = obs;
        st.candidate_run = 1;
        st.candidate_start_frame = frame_index;
        st.candidate_start_timestep = frame.timestep;
        continue;
      }

      if (obs == st.candidate_site) {
        ++st.candidate_run;
      } else {
        st.candidate_site = obs;
        st.candidate_run = 1;
        st.candidate_start_frame = frame_index;
        st.candidate_start_timestep = frame.timestep;
      }

      if (st.candidate_site != st.committed_site && st.candidate_run >= opt_.min_persistence_frames) {
        commit_transition_(st, st.candidate_site, st.candidate_start_frame, st.candidate_start_timestep);
      }
    }

    ++n_frames_used_;
    if (n_frames_used_ == 1) {
      first_timestep_used_ = frame.timestep;
      first_frame_index_used_ = frame_index;
    }
    last_timestep_used_ = frame.timestep;
    last_frame_index_used_ = frame_index;
  }

  void flush_partial() override {
    write_outputs_();
  }

  void finalize() override {
    write_outputs_();
  }

private:
  struct ParticleState {
    bool initialized = false;
    std::int64_t committed_site = 0;
    std::size_t entry_frame = 0;
    std::int64_t entry_timestep = 0;

    std::int64_t candidate_site = 0;
    std::size_t candidate_run = 0;
    std::size_t candidate_start_frame = 0;
    std::int64_t candidate_start_timestep = 0;
  };

  void validate_() const {
    if (sel_.idx.empty()) throw std::runtime_error("site_hop_stats: selection is empty");
    if (opt_.range.frame_start < 0) throw std::runtime_error("site_hop_stats: frame_start must be >= 0");
    if (opt_.range.frame_end >= 0 && opt_.range.frame_end < opt_.range.frame_start) {
      throw std::runtime_error("site_hop_stats: frame_end must be -1 or >= frame_start");
    }
    if (opt_.min_persistence_frames == 0) {
      throw std::runtime_error("site_hop_stats: min_persistence_frames must be >= 1");
    }
  }

  bool frame_in_range_(std::size_t frame_index) const {
    const std::int64_t fi = static_cast<std::int64_t>(frame_index);
    if (fi < opt_.range.frame_start) return false;
    if (opt_.range.frame_end >= 0 && fi > opt_.range.frame_end) return false;
    return true;
  }

  bool valid_site_(std::int64_t s) const {
    return !(opt_.ignore_unassigned && s == opt_.unassigned_site);
  }

  void commit_transition_(ParticleState& st,
                          std::int64_t new_site,
                          std::size_t new_entry_frame,
                          std::int64_t new_entry_timestep) {
    const std::size_t dwell_frames = (new_entry_frame >= st.entry_frame) ? (new_entry_frame - st.entry_frame) : 0;
    const double dwell_time = (new_entry_timestep >= st.entry_timestep)
        ? (opt_.dt * static_cast<double>(new_entry_timestep - st.entry_timestep))
        : 0.0;

    if (valid_site_(st.committed_site)) {
      if (dwell_frames > 0) {
        dwell_hist_[dwell_frames] += 1;
        dwell_time_sum_ += dwell_time;
        dwell_frame_sum_ += static_cast<double>(dwell_frames);
        ++n_complete_dwells_;
      }
    }

    if (valid_site_(st.committed_site) && valid_site_(new_site)) {
      ++transition_counts_[{st.committed_site, new_site}];
      ++n_hops_total_;
    }

    st.committed_site = new_site;
    st.entry_frame = new_entry_frame;
    st.entry_timestep = new_entry_timestep;
  }

  void write_outputs_() const {
    write_summary_();
    write_transitions_();
    write_dwell_();
  }

  void write_summary_() const {
    util::atomic_write_text(summary_path_, [&](std::ostream& ofs) {
      ofs << std::setprecision(17);
      ofs << "# PILOTS: site hop summary\n";
      ofs << "# selection: " << sel_.name << "\n";
      ofs << "# site_field: " << opt_.site_field
          << ", unassigned_site: " << opt_.unassigned_site
          << ", min_persistence_frames: " << opt_.min_persistence_frames
          << ", ignore_unassigned: " << (opt_.ignore_unassigned ? 1 : 0) << "\n";
      ofs << "# columns: n_particles  n_frames  n_hops_total  hop_rate_per_particle_per_frame  hop_rate_per_particle_per_time  mean_dwell_frames  mean_dwell_time\n";
      const double obs_frames = (n_frames_used_ > 0) ? static_cast<double>(n_frames_used_) : 0.0;
      const double obs_time = (n_frames_used_ > 1)
          ? (opt_.dt * static_cast<double>(last_timestep_used_ - first_timestep_used_))
          : 0.0;
      const double n_particles = static_cast<double>(sel_.idx.size());
      const double hop_rate_frame = (obs_frames > 0.0 && n_particles > 0.0)
          ? (static_cast<double>(n_hops_total_) / (n_particles * obs_frames))
          : 0.0;
      const double hop_rate_time = (obs_time > 0.0 && n_particles > 0.0)
          ? (static_cast<double>(n_hops_total_) / (n_particles * obs_time))
          : 0.0;
      const double mean_dwell_frames = (n_complete_dwells_ > 0)
          ? (dwell_frame_sum_ / static_cast<double>(n_complete_dwells_))
          : 0.0;
      const double mean_dwell_time = (n_complete_dwells_ > 0)
          ? (dwell_time_sum_ / static_cast<double>(n_complete_dwells_))
          : 0.0;
      ofs << sel_.idx.size() << ' '
          << n_frames_used_ << ' '
          << n_hops_total_ << ' '
          << hop_rate_frame << ' '
          << hop_rate_time << ' '
          << mean_dwell_frames << ' '
          << mean_dwell_time << '\n';
    });
  }

  void write_transitions_() const {
    util::atomic_write_text(transitions_path_, [&](std::ostream& ofs) {
      ofs << std::setprecision(17);
      ofs << "# PILOTS: site hop transition counts\n";
      ofs << "# columns: site_from  site_to  count  probability\n";
      for (const auto& kv : transition_counts_) {
        const double prob = (n_hops_total_ > 0)
            ? (static_cast<double>(kv.second) / static_cast<double>(n_hops_total_))
            : 0.0;
        ofs << kv.first.first << ' ' << kv.first.second << ' ' << kv.second << ' ' << prob << '\n';
      }
    });
  }

  void write_dwell_() const {
    util::atomic_write_text(dwell_path_, [&](std::ostream& ofs) {
      ofs << std::setprecision(17);
      ofs << "# PILOTS: site dwell histogram (completed dwells only)\n";
      ofs << "# columns: dwell_frames  dwell_time  count  probability\n";
      for (const auto& kv : dwell_hist_) {
        const double dwell_time = opt_.dt * static_cast<double>(kv.first);
        const double prob = (n_complete_dwells_ > 0)
            ? (static_cast<double>(kv.second) / static_cast<double>(n_complete_dwells_))
            : 0.0;
        ofs << kv.first << ' ' << dwell_time << ' ' << kv.second << ' ' << prob << '\n';
      }
    });
  }

  std::string instance_name_;
  std::string summary_path_;
  std::string transitions_path_;
  std::string dwell_path_;
  SiteHopOptions opt_;
  std::string sel_name_owned_;
  SelectionView sel_{};
  std::vector<ParticleState> states_;
  std::map<std::pair<std::int64_t, std::int64_t>, std::size_t> transition_counts_;
  std::map<std::size_t, std::size_t> dwell_hist_;
  std::size_t n_hops_total_ = 0;
  double dwell_frame_sum_ = 0.0;
  double dwell_time_sum_ = 0.0;
  std::size_t n_complete_dwells_ = 0;
  std::size_t n_frames_used_ = 0;
  std::int64_t first_timestep_used_ = 0;
  std::int64_t last_timestep_used_ = 0;
  std::size_t first_frame_index_used_ = 0;
  std::size_t last_frame_index_used_ = 0;
};

MeasureCapabilities current_acf_caps(const IniConfig& cfg,
                                     const std::string& section,
                                     const std::string& instance,
                                     const MeasureBuildEnv& env) {
  (void)instance;
  (void)env;
  MeasureCapabilities caps;
  append_static_selection_caps(cfg, section, caps);
  caps.requires_dfields = {
      cfg.get_string(section, "q_field", std::optional<std::string>("q")),
      cfg.get_string(section, "vector_x_field", std::optional<std::string>("vx")),
      cfg.get_string(section, "vector_y_field", std::optional<std::string>("vy")),
      cfg.get_string(section, "vector_z_field", std::optional<std::string>("vz"))};
  return caps;
}

std::unique_ptr<IMeasure> current_acf_create(const IniConfig& cfg,
                                             const std::string& section,
                                             const std::string& instance,
                                             const MeasureBuildEnv& env,
                                             const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("current_acf factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "current_acf");

  CurrentACFOptions opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.diag_mask = parse_diag_mask(cfg.get_string(section, "components", std::optional<std::string>("xxyyzz")));
  opt.q_field = cfg.get_string(section, "q_field", std::optional<std::string>("q"));
  opt.vector_x_field = cfg.get_string(section, "vector_x_field", std::optional<std::string>("vx"));
  opt.vector_y_field = cfg.get_string(section, "vector_y_field", std::optional<std::string>("vy"));
  opt.vector_z_field = cfg.get_string(section, "vector_z_field", std::optional<std::string>("vz"));
  opt.kbt = cfg.get_double(section, "kbt");
  opt.corr = parse_correlator_spec(cfg, section, env.dt);
  opt.range.frame_end = resolve_exact_frame_end(opt.corr, env.follow, opt.range.frame_start, opt.range.frame_end, env.input_path);

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("current_acf.dat"))).lexically_normal();

  return std::make_unique<CurrentACFMeasure>(instance, out.string(), sel, opt);
}

MeasureCapabilities conductivity_einstein_caps(const IniConfig& cfg,
                                               const std::string& section,
                                               const std::string& instance,
                                               const MeasureBuildEnv& env) {
  (void)instance;
  (void)env;
  MeasureCapabilities caps;
  append_static_selection_caps(cfg, section, caps);
  caps.requires_dfields = {"xu", "yu", "zu", cfg.get_string(section, "q_field", std::optional<std::string>("q"))};
  return caps;
}

std::unique_ptr<IMeasure> conductivity_einstein_create(const IniConfig& cfg,
                                                       const std::string& section,
                                                       const std::string& instance,
                                                       const MeasureBuildEnv& env,
                                                       const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("conductivity_einstein factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "conductivity_einstein");

  ChargeEinsteinOptions opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.diag_mask = parse_diag_mask(cfg.get_string(section, "components", std::optional<std::string>("xxyyzz")));
  opt.q_field = cfg.get_string(section, "q_field", std::optional<std::string>("q"));
  opt.kbt = cfg.get_double(section, "kbt");
  opt.corr = parse_correlator_spec(cfg, section, env.dt);
  opt.range.frame_end = resolve_exact_frame_end(opt.corr, env.follow, opt.range.frame_start, opt.range.frame_end, env.input_path);

  const auto q0 = frame0.require_dfield(opt.q_field);
  const std::vector<double> q_sel = selected_field_copy(sel, q0);

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("conductivity_einstein.dat"))).lexically_normal();

  return std::make_unique<ChargeEinsteinMeasure>(instance, out.string(), sel, opt, q_sel);
}

MeasureCapabilities site_hop_caps(const IniConfig& cfg,
                                  const std::string& section,
                                  const std::string& instance,
                                  const MeasureBuildEnv& env) {
  (void)instance;
  (void)env;
  MeasureCapabilities caps;
  append_static_selection_caps(cfg, section, caps);
  append_integer_like_field_cap(caps, cfg.get_string(section, "site_field", std::optional<std::string>("site")));
  return caps;
}

std::unique_ptr<IMeasure> site_hop_create(const IniConfig& cfg,
                                          const std::string& section,
                                          const std::string& instance,
                                          const MeasureBuildEnv& env,
                                          const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("site_hop_stats factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "site_hop_stats");

  SiteHopOptions opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.site_field = cfg.get_string(section, "site_field", std::optional<std::string>("site"));
  opt.unassigned_site = cfg.get_int64(section, "unassigned_site", std::optional<std::int64_t>(-1));
  opt.min_persistence_frames = cfg.get_size(section, "min_persistence_frames", std::optional<std::size_t>(1));
  opt.ignore_unassigned = cfg.get_bool(section, "ignore_unassigned", std::optional<bool>(true));
  opt.dt = env.dt;

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path summary = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("site_hop_summary.dat"))).lexically_normal();
  const fs::path transitions = (output_dir / cfg.get_string(section, "transitions_output", std::optional<std::string>("site_hop_transitions.dat"))).lexically_normal();
  const fs::path dwell = (output_dir / cfg.get_string(section, "dwell_output", std::optional<std::string>("site_hop_dwell.dat"))).lexically_normal();

  return std::make_unique<SiteHopStatsMeasure>(instance, summary.string(), transitions.string(), dwell.string(), sel, opt);
}

static MeasureRegistrar g_register_current_acf("current_acf", &current_acf_caps, &current_acf_create);
static MeasureRegistrar g_register_conductivity_gk_alias("conductivity_gk", &current_acf_caps, &current_acf_create);
static MeasureRegistrar g_register_conductivity_einstein("conductivity_einstein", &conductivity_einstein_caps, &conductivity_einstein_create);
static MeasureRegistrar g_register_ionic_conductivity_einstein_alias("ionic_conductivity_einstein", &conductivity_einstein_caps, &conductivity_einstein_create);
static MeasureRegistrar g_register_charge_msd_alias("charge_msd", &conductivity_einstein_caps, &conductivity_einstein_create);
static MeasureRegistrar g_register_site_hop_stats("site_hop_stats", &site_hop_caps, &site_hop_create);
static MeasureRegistrar g_register_hop_stats_alias("hop_stats", &site_hop_caps, &site_hop_create);

} // namespace
} // namespace pilots
