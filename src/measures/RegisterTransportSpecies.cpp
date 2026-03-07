#include <algorithm>
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
#include <unordered_map>
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

using measure_ext::Vec3;
using measure_ext::count_diag_dims;
using measure_ext::dstr;
using measure_ext::get_static_combined_view;
using measure_ext::get_static_group_view;
using measure_ext::integer_like_field_to_i64;
using measure_ext::lag_axis_name;
using measure_ext::resolve_exact_frame_end;
using measure_ext::resolve_path;
using measure_ext::x_unit_for_axis;

struct RangeOptions {
  std::int64_t frame_start = 0;
  std::int64_t frame_end = -1;
  bool dry_run = false;
};

struct LagOptions {
  std::int64_t max_lag_frames = -1;
  std::int64_t lag_stride = 1;
  std::int64_t origin_stride = 1;
};

inline bool frame_in_range(std::size_t frame_index, const RangeOptions& opt) {
  const std::int64_t fi = static_cast<std::int64_t>(frame_index);
  if (fi < opt.frame_start) return false;
  if (opt.frame_end >= 0 && fi > opt.frame_end) return false;
  return true;
}

inline LagOptions lag_options_from_config(const IniConfig& cfg,
                                          const std::string& section) {
  LagOptions opt;
  opt.max_lag_frames = cfg.get_int64(section, "max_lag_frames", std::optional<std::int64_t>(-1));
  opt.lag_stride = cfg.get_int64(section, "lag_stride", std::optional<std::int64_t>(1));
  opt.origin_stride = cfg.get_int64(section, "origin_stride", std::optional<std::int64_t>(1));
  if (opt.max_lag_frames == 0 || opt.max_lag_frames < -1) {
    throw std::runtime_error(section + ": max_lag_frames must be -1 or >= 1");
  }
  if (opt.lag_stride <= 0) throw std::runtime_error(section + ": lag_stride must be >= 1");
  if (opt.origin_stride <= 0) throw std::runtime_error(section + ": origin_stride must be >= 1");
  return opt;
}

inline std::vector<std::size_t> build_lag_list(std::size_t nframes,
                                               const LagOptions& opt,
                                               bool include_zero = true) {
  std::vector<std::size_t> out;
  if (nframes == 0) return out;
  const std::size_t max_lag = (opt.max_lag_frames < 0)
      ? (nframes - 1)
      : std::min<std::size_t>(static_cast<std::size_t>(opt.max_lag_frames), nframes - 1);
  std::size_t start = include_zero ? 0 : 1;
  for (std::size_t lag = start; lag <= max_lag; lag += static_cast<std::size_t>(opt.lag_stride)) {
    out.push_back(lag);
  }
  if (include_zero && out.empty()) out.push_back(0);
  return out;
}

inline double mean_dtimestep_for_lag(const std::vector<std::int64_t>& timesteps,
                                     std::size_t lag,
                                     std::size_t origin_stride) {
  if (timesteps.empty() || lag >= timesteps.size()) return 0.0;
  double sum = 0.0;
  std::size_t n = 0;
  for (std::size_t o = 0; o + lag < timesteps.size(); o += origin_stride) {
    sum += static_cast<double>(timesteps[o + lag] - timesteps[o]);
    ++n;
  }
  return (n > 0) ? (sum / static_cast<double>(n)) : 0.0;
}

struct SpeciesGroups {
  std::vector<std::int64_t> species_id;
  std::vector<std::vector<std::size_t>> selpos_by_species;
};

inline SpeciesGroups build_species_groups(const SelectionView& sel,
                                          std::span<const std::int64_t> species_per_atom) {
  std::unordered_map<std::int64_t, std::vector<std::size_t>> tmp;
  tmp.reserve(sel.idx.size() / 4 + 8);
  for (std::size_t p = 0; p < sel.idx.size(); ++p) {
    const std::size_t atom = sel.idx[p];
    if (atom >= species_per_atom.size()) {
      throw std::runtime_error("species field vector smaller than selected atom index");
    }
    tmp[species_per_atom[atom]].push_back(p);
  }
  std::vector<std::int64_t> ids;
  ids.reserve(tmp.size());
  for (const auto& kv : tmp) ids.push_back(kv.first);
  std::sort(ids.begin(), ids.end());

  SpeciesGroups out;
  out.species_id = ids;
  out.selpos_by_species.reserve(ids.size());
  for (const auto sid : ids) out.selpos_by_species.push_back(std::move(tmp[sid]));
  return out;
}

inline std::vector<double> selected_field_copy(const SelectionView& sel,
                                               std::span<const double> field) {
  std::vector<double> out(sel.idx.size(), 0.0);
#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
  for (std::size_t p = 0; p < sel.idx.size(); ++p) out[p] = field[sel.idx[p]];
  return out;
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
      sxx += dx * dx;
      syy += dy * dy;
      szz += dz * dz;
      sxy += dx * dy;
      sxz += dx * dz;
      syz += dy * dz;
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

struct VACFPairOpT6 {
  Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const {
    const std::size_t n = cur.xx.size();
    if (org.xx.size() != n || cur.yy.size() != n || org.yy.size() != n ||
        cur.zz.size() != n || org.zz.size() != n) {
      throw std::runtime_error("VACFPairOpT6: slot sizes mismatch");
    }

    double sxx = 0.0, syy = 0.0, szz = 0.0, sxy = 0.0, sxz = 0.0, syz = 0.0;
#if PILOTS_HAS_OPENMP
#pragma omp parallel for reduction(+:sxx,syy,szz,sxy,sxz,syz)
#endif
    for (std::size_t p = 0; p < n; ++p) {
      const double cx = cur.xx[p];
      const double cy = cur.yy[p];
      const double cz = cur.zz[p];
      const double ox = org.xx[p];
      const double oy = org.yy[p];
      const double oz = org.zz[p];
      sxx += cx * ox;
      syy += cy * oy;
      szz += cz * oz;
      sxy += cx * oy;
      sxz += cx * oz;
      syz += cy * oz;
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

enum class SpeciesTensorMode {
  MSD,
  VACF,
};

inline const char* species_tensor_mode_name(SpeciesTensorMode mode) {
  switch (mode) {
    case SpeciesTensorMode::MSD: return "species_msd";
    case SpeciesTensorMode::VACF: return "vacf";
  }
  return "species_msd";
}

struct SpeciesTensorOptions {
  RangeOptions range;
  CorrelatorSpec corr;
  SpeciesTensorMode mode = SpeciesTensorMode::MSD;
  bool remove_drift = false;
  std::string vector_x_field = "vx";
  std::string vector_y_field = "vy";
  std::string vector_z_field = "vz";
};

class SpeciesTensorCorrMeasure final : public IMeasure {
public:
  SpeciesTensorCorrMeasure(std::string type_name,
                           std::string instance_name,
                           std::string output_path,
                           SelectionView sel,
                           SelectionView drift_sel,
                           SpeciesTensorOptions opt,
                           SpeciesGroups species)
      : type_name_(std::move(type_name)),
        instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)),
        species_(std::move(species)) {
    sel_name_owned_ = std::string(sel.name);
    drift_name_owned_ = std::string(drift_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    drift_sel_ = SelectionView{drift_name_owned_, drift_sel.idx};
    validate_();

    for (const auto& idx : species_.selpos_by_species) {
      const std::size_t n = idx.size();
      if (opt_.mode == SpeciesTensorMode::MSD) {
        corr_.push_back(make_correlator_t6_runtime_auto(n, exact_window_frames_(), opt_.corr, CartesianMSDPairOpT6{}));
      } else {
        corr_.push_back(make_correlator_t6_runtime_auto(n, exact_window_frames_(), opt_.corr, VACFPairOpT6{}));
      }
      scratch_xx_.push_back(std::vector<double>(n, 0.0));
      scratch_yy_.push_back(std::vector<double>(n, 0.0));
      scratch_zz_.push_back(std::vector<double>(n, 0.0));
      scratch_zero_.push_back(std::vector<double>(n, 0.0));
    }
  }

  std::string type() const override { return type_name_; }
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
    od.columns = {"lag", "time", "species_id",
                  "xx", "yy", "zz", "xy", "xz", "yz", "trace",
                  "aux_xx", "aux_yy", "aux_zz", "aux_trace",
                  "count_pairs", "n_blocks", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));

    md.params["mode"] = type_name_;
    md.params["frame_start"] = std::to_string(opt_.range.frame_start);
    md.params["frame_end"] = std::to_string(opt_.range.frame_end);
    md.params["correlator"] = opt_.corr.type;
    md.params["lag_axis"] = lag_axis_name(opt_.corr.axis);
    md.params["remove_drift"] = opt_.remove_drift ? "true" : "false";
    if (opt_.mode == SpeciesTensorMode::VACF) {
      md.params["vector_x_field"] = opt_.vector_x_field;
      md.params["vector_y_field"] = opt_.vector_y_field;
      md.params["vector_z_field"] = opt_.vector_z_field;
    }
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

    double comx = 0.0, comy = 0.0, comz = 0.0;
    if (opt_.remove_drift) {
      const auto xu = frame.require_dfield("xu");
      const auto yu = frame.require_dfield("yu");
      const auto zu = frame.require_dfield("zu");
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

    std::span<const double> fx;
    std::span<const double> fy;
    std::span<const double> fz;
    if (opt_.mode == SpeciesTensorMode::MSD) {
      fx = frame.require_dfield("xu");
      fy = frame.require_dfield("yu");
      fz = frame.require_dfield("zu");
    } else {
      fx = frame.require_dfield(opt_.vector_x_field);
      fy = frame.require_dfield(opt_.vector_y_field);
      fz = frame.require_dfield(opt_.vector_z_field);
    }

    for (std::size_t s = 0; s < species_.species_id.size(); ++s) {
      const auto& idx = species_.selpos_by_species[s];
#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
      for (std::size_t k = 0; k < idx.size(); ++k) {
        const std::size_t p = idx[k];
        const std::size_t atom = sel_.idx[p];
        scratch_xx_[s][k] = fx[atom] - comx;
        scratch_yy_[s][k] = fy[atom] - comy;
        scratch_zz_[s][k] = fz[atom] - comz;
      }
      corr_[s]->push(frame.timestep,
                     std::span<const double>(scratch_xx_[s].data(), scratch_xx_[s].size()),
                     std::span<const double>(scratch_yy_[s].data(), scratch_yy_[s].size()),
                     std::span<const double>(scratch_zz_[s].data(), scratch_zz_[s].size()),
                     std::span<const double>(scratch_zero_[s].data(), scratch_zero_[s].size()),
                     std::span<const double>(scratch_zero_[s].data(), scratch_zero_[s].size()),
                     std::span<const double>(scratch_zero_[s].data(), scratch_zero_[s].size()));
    }
  }

  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: species-resolved " << type_name_ << "\n";
      ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
      ofs << "# frame_range: [" << opt_.range.frame_start << ", " << opt_.range.frame_end << "]\n";
      ofs << "# correlator: " << opt_.corr.type << ", lag_axis: " << lag_axis_name(opt_.corr.axis) << "\n";
      ofs << "# remove_drift: " << (opt_.remove_drift ? "true" : "false") << ", drift_group: " << drift_sel_.name << "\n";
      if (opt_.mode == SpeciesTensorMode::VACF) {
        ofs << "# velocity_fields: " << opt_.vector_x_field << ", " << opt_.vector_y_field << ", " << opt_.vector_z_field << "\n";
      }
      ofs << "# columns: lag time species_id xx yy zz xy xz yz trace aux_xx aux_yy aux_zz aux_trace count_pairs n_blocks mean_dtimestep\n";
      for (std::size_t s = 0; s < species_.species_id.size(); ++s) {
        const CorrelationSeriesT6 series = corr_[s]->snapshot();
        for (std::size_t i = 0; i < series.lag.size(); ++i) {
          const double xx = series.value_xx[i];
          const double yy = series.value_yy[i];
          const double zz = series.value_zz[i];
          const double xy = series.value_xy[i];
          const double xz = series.value_xz[i];
          const double yz = series.value_yz[i];
          const double trace = xx + yy + zz;
          double aux_xx = 0.0, aux_yy = 0.0, aux_zz = 0.0, aux_trace = 0.0;
          if (opt_.mode == SpeciesTensorMode::MSD) {
            const double t = series.time[i];
            if (t > 0.0) {
              aux_xx = 0.5 * xx / t;
              aux_yy = 0.5 * yy / t;
              aux_zz = 0.5 * zz / t;
              aux_trace = trace / (6.0 * t);
            }
          } else {
            // running diffusion tensor from VACF integral is intentionally not emitted here;
            // the aux columns are reserved to keep a common schema across aliases.
          }
          ofs << std::setprecision(17) << series.lag[i] << ' '
              << std::setprecision(17) << series.time[i] << ' '
              << species_.species_id[s] << ' '
              << std::setprecision(17) << xx << ' '
              << std::setprecision(17) << yy << ' '
              << std::setprecision(17) << zz << ' '
              << std::setprecision(17) << xy << ' '
              << std::setprecision(17) << xz << ' '
              << std::setprecision(17) << yz << ' '
              << std::setprecision(17) << trace << ' '
              << std::setprecision(17) << aux_xx << ' '
              << std::setprecision(17) << aux_yy << ' '
              << std::setprecision(17) << aux_zz << ' '
              << std::setprecision(17) << aux_trace << ' '
              << series.count_pairs[i] << ' '
              << series.n_blocks[i] << ' '
              << std::setprecision(17) << series.mean_dtimestep[i] << '\n';
        }
      }
    });
  }

  void finalize() override {
    if (opt_.range.dry_run) return;
    flush_partial();
  }

private:
  void validate_() const {
    if (sel_.idx.empty()) throw std::runtime_error(type_name_ + ": selection is empty");
    if (species_.species_id.empty()) throw std::runtime_error(type_name_ + ": no species found in selection");
    if (opt_.range.frame_start < 0) throw std::runtime_error(type_name_ + ": frame_start must be >= 0");
    if (opt_.range.frame_end >= 0 && opt_.range.frame_end < opt_.range.frame_start) {
      throw std::runtime_error(type_name_ + ": frame_end must be -1 or >= frame_start");
    }
  }

  std::size_t exact_window_frames_() const {
    if (opt_.corr.type != "exact") return 0;
    if (opt_.range.frame_end < 0) {
      throw std::runtime_error(type_name_ + ": exact correlator requires finite frame_end (factory must resolve this)");
    }
    return static_cast<std::size_t>(opt_.range.frame_end - opt_.range.frame_start + 1);
  }

  std::string type_name_;
  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  std::string drift_name_owned_;
  SelectionView sel_;
  SelectionView drift_sel_;
  SpeciesTensorOptions opt_;
  SpeciesGroups species_;
  bool started_ = false;

  std::vector<std::unique_ptr<ICorrelatorT6>> corr_;
  std::vector<std::vector<double>> scratch_xx_, scratch_yy_, scratch_zz_, scratch_zero_;
};

class SpeciesCurrentCorrMeasure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    LagOptions lag;
    std::string q_field = "q";
    std::string vx_field = "vx";
    std::string vy_field = "vy";
    std::string vz_field = "vz";
  };

  SpeciesCurrentCorrMeasure(std::string instance_name,
                            std::string output_path,
                            SelectionView sel,
                            SpeciesGroups species,
                            Options opt)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        species_(std::move(species)),
        opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    if (sel_.idx.empty()) throw std::runtime_error("species_current_corr: selection is empty");
    if (species_.species_id.empty()) throw std::runtime_error("species_current_corr: no species in selection");
  }

  std::string type() const override { return "species_current_corr"; }
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
    od.x_axis = "lag";
    od.x_unit = "frames";
    od.columns = {"lag", "time", "species_i", "species_j",
                  "c_xx", "c_yy", "c_zz", "c_xy", "c_xz", "c_yz", "c_trace",
                  "count_pairs", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));

    md.params["q_field"] = opt_.q_field;
    md.params["vx_field"] = opt_.vx_field;
    md.params["vy_field"] = opt_.vy_field;
    md.params["vz_field"] = opt_.vz_field;
    return md;
  }

  void on_start(const Frame& first_frame) override {
    const auto q = first_frame.require_dfield(opt_.q_field);
    q_sel_.resize(sel_.idx.size(), 0.0);
    for (std::size_t p = 0; p < sel_.idx.size(); ++p) q_sel_[p] = q[sel_.idx[p]];
    started_ = true;
    if (!opt_.range.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    const auto vx = frame.require_dfield(opt_.vx_field);
    const auto vy = frame.require_dfield(opt_.vy_field);
    const auto vz = frame.require_dfield(opt_.vz_field);

    std::vector<Vec3> js(species_.species_id.size(), Vec3{});
    for (std::size_t s = 0; s < species_.species_id.size(); ++s) {
      double jx = 0.0, jy = 0.0, jz = 0.0;
      for (const std::size_t p : species_.selpos_by_species[s]) {
        const std::size_t atom = sel_.idx[p];
        const double qi = q_sel_[p];
        jx += qi * vx[atom];
        jy += qi * vy[atom];
        jz += qi * vz[atom];
      }
      js[s] = Vec3{jx, jy, jz};
    }
    series_.push_back(std::move(js));
    timesteps_.push_back(frame.timestep);
  }

  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    write_output_();
  }

  void finalize() override {
    if (opt_.range.dry_run) return;
    write_output_();
  }

private:
  void write_output_() {
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: species-resolved current cross-correlation matrix\n";
      ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
      ofs << "# columns: lag time species_i species_j c_xx c_yy c_zz c_xy c_xz c_yz c_trace count_pairs mean_dtimestep\n";
      const auto lags = build_lag_list(series_.size(), opt_.lag, true);
      for (const std::size_t lag : lags) {
        const double time = mean_dtimestep_for_lag(timesteps_, lag, static_cast<std::size_t>(opt_.lag.origin_stride));
        for (std::size_t a = 0; a < species_.species_id.size(); ++a) {
          for (std::size_t b = 0; b < species_.species_id.size(); ++b) {
            double sxx = 0.0, syy = 0.0, szz = 0.0, sxy = 0.0, sxz = 0.0, syz = 0.0;
            std::size_t n = 0;
            for (std::size_t o = 0; o + lag < series_.size(); o += static_cast<std::size_t>(opt_.lag.origin_stride)) {
              const Vec3& ja = series_[o + lag][a];
              const Vec3& jb = series_[o][b];
              sxx += ja.x * jb.x;
              syy += ja.y * jb.y;
              szz += ja.z * jb.z;
              sxy += ja.x * jb.y;
              sxz += ja.x * jb.z;
              syz += ja.y * jb.z;
              ++n;
            }
            const double inv = (n > 0) ? (1.0 / static_cast<double>(n)) : 0.0;
            const double cxx = sxx * inv;
            const double cyy = syy * inv;
            const double czz = szz * inv;
            const double cxy = sxy * inv;
            const double cxz = sxz * inv;
            const double cyz = syz * inv;
            ofs << lag << ' '
                << std::setprecision(17) << time << ' '
                << species_.species_id[a] << ' '
                << species_.species_id[b] << ' '
                << std::setprecision(17) << cxx << ' '
                << std::setprecision(17) << cyy << ' '
                << std::setprecision(17) << czz << ' '
                << std::setprecision(17) << cxy << ' '
                << std::setprecision(17) << cxz << ' '
                << std::setprecision(17) << cyz << ' '
                << std::setprecision(17) << (cxx + cyy + czz) << ' '
                << n << ' '
                << std::setprecision(17) << time << '\n';
          }
        }
      }
    });
  }

  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  SelectionView sel_;
  SpeciesGroups species_;
  Options opt_;
  bool started_ = false;
  std::vector<double> q_sel_;
  std::vector<std::vector<Vec3>> series_;
  std::vector<std::int64_t> timesteps_;
};

class OnsagerTransportMeasure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    LagOptions lag;
    std::string q_field = "q";
    bool remove_drift = false;
    double kbt = 1.0;
  };

  OnsagerTransportMeasure(std::string type_name,
                          std::string instance_name,
                          std::string matrix_output_path,
                          std::string transference_output_path,
                          SelectionView sel,
                          SelectionView drift_sel,
                          SpeciesGroups species,
                          Options opt)
      : type_name_(std::move(type_name)),
        instance_name_(std::move(instance_name)),
        matrix_output_path_(std::move(matrix_output_path)),
        transference_output_path_(std::move(transference_output_path)),
        species_(std::move(species)),
        opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    drift_name_owned_ = std::string(drift_sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    drift_sel_ = SelectionView{drift_name_owned_, drift_sel.idx};
    if (sel_.idx.empty()) throw std::runtime_error(type_name_ + ": selection is empty");
    if (!(opt_.kbt > 0.0)) throw std::runtime_error(type_name_ + ": kbt must be > 0");
  }

  std::string type() const override { return type_name_; }
  std::string instance_name() const override { return instance_name_; }

  output::MeasureDescriptor describe() const override {
    output::MeasureDescriptor md;
    md.instance = instance_name_;
    md.type = type();
    md.selection = std::string(sel_.name);
    md.n_selected = sel_.idx.size();

    output::OutputFileDescriptor od1;
    od1.path = matrix_output_path_;
    od1.format = "text";
    od1.x_axis = "lag";
    od1.x_unit = "frames";
    od1.columns = {"lag", "time", "species_i", "species_j",
                   "l_xx", "l_yy", "l_zz", "l_xy", "l_xz", "l_yz", "l_trace",
                   "count_pairs", "mean_dtimestep"};
    md.outputs.push_back(std::move(od1));

    output::OutputFileDescriptor od2;
    od2.path = transference_output_path_;
    od2.format = "text";
    od2.x_axis = "lag";
    od2.x_unit = "frames";
    od2.columns = {"lag", "time", "species_id", "t_plus_like", "row_sum_trace", "total_trace", "count_pairs", "mean_dtimestep"};
    md.outputs.push_back(std::move(od2));

    md.params["q_field"] = opt_.q_field;
    md.params["kbt"] = dstr(opt_.kbt);
    md.params["remove_drift"] = opt_.remove_drift ? "true" : "false";
    return md;
  }

  void on_start(const Frame& first_frame) override {
    const auto q = first_frame.require_dfield(opt_.q_field);
    q_sel_.resize(sel_.idx.size(), 0.0);
    for (std::size_t p = 0; p < sel_.idx.size(); ++p) q_sel_[p] = q[sel_.idx[p]];
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
      const std::size_t nd = drift_sel_.idx.size();
      for (std::size_t k = 0; k < nd; ++k) {
        const std::size_t i = drift_sel_.idx[k];
        sx += xu[i]; sy += yu[i]; sz += zu[i];
      }
      const double inv = 1.0 / static_cast<double>(nd);
      comx = sx * inv; comy = sy * inv; comz = sz * inv;
    }

    x_series_.push_back(std::vector<double>(sel_.idx.size(), 0.0));
    y_series_.push_back(std::vector<double>(sel_.idx.size(), 0.0));
    z_series_.push_back(std::vector<double>(sel_.idx.size(), 0.0));
    auto& xs = x_series_.back();
    auto& ys = y_series_.back();
    auto& zs = z_series_.back();
    for (std::size_t p = 0; p < sel_.idx.size(); ++p) {
      const std::size_t atom = sel_.idx[p];
      xs[p] = xu[atom] - comx;
      ys[p] = yu[atom] - comy;
      zs[p] = zu[atom] - comz;
    }
    timesteps_.push_back(frame.timestep);
    volume_sum_ += frame.box.lx() * frame.box.ly() * frame.box.lz();
  }

  void flush_partial() override {
    if (opt_.range.dry_run || !started_) return;
    write_outputs_();
  }

  void finalize() override {
    if (opt_.range.dry_run) return;
    write_outputs_();
  }

private:
  void write_outputs_() {
    const double mean_volume = (!timesteps_.empty()) ? (volume_sum_ / static_cast<double>(timesteps_.size())) : 0.0;
    const auto lags = build_lag_list(timesteps_.size(), opt_.lag, false);

    util::atomic_write_text(matrix_output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: Einstein-form species Onsager transport matrix\n";
      ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
      ofs << "# remove_drift: " << (opt_.remove_drift ? "true" : "false") << ", drift_group: " << drift_sel_.name << "\n";
      ofs << "# kBT: " << std::setprecision(17) << opt_.kbt << ", mean_volume: " << std::setprecision(17) << mean_volume << "\n";
      ofs << "# columns: lag time species_i species_j l_xx l_yy l_zz l_xy l_xz l_yz l_trace count_pairs mean_dtimestep\n";

      for (const std::size_t lag : lags) {
        const double dt = mean_dtimestep_for_lag(timesteps_, lag, static_cast<std::size_t>(opt_.lag.origin_stride));
        const double denom_tensor = (dt > 0.0 && mean_volume > 0.0) ? (2.0 * mean_volume * opt_.kbt * dt) : std::numeric_limits<double>::quiet_NaN();
        for (std::size_t a = 0; a < species_.species_id.size(); ++a) {
          for (std::size_t b = 0; b < species_.species_id.size(); ++b) {
            double sxx = 0.0, syy = 0.0, szz = 0.0, sxy = 0.0, sxz = 0.0, syz = 0.0;
            std::size_t n = 0;
            for (std::size_t o = 0; o + lag < x_series_.size(); o += static_cast<std::size_t>(opt_.lag.origin_stride)) {
              double dmax = 0.0, dmay = 0.0, dmaz = 0.0;
              double dmbx = 0.0, dmby = 0.0, dmbz = 0.0;
              for (const std::size_t p : species_.selpos_by_species[a]) {
                const double qi = q_sel_[p];
                dmax += qi * (x_series_[o + lag][p] - x_series_[o][p]);
                dmay += qi * (y_series_[o + lag][p] - y_series_[o][p]);
                dmaz += qi * (z_series_[o + lag][p] - z_series_[o][p]);
              }
              for (const std::size_t p : species_.selpos_by_species[b]) {
                const double qi = q_sel_[p];
                dmbx += qi * (x_series_[o + lag][p] - x_series_[o][p]);
                dmby += qi * (y_series_[o + lag][p] - y_series_[o][p]);
                dmbz += qi * (z_series_[o + lag][p] - z_series_[o][p]);
              }
              sxx += dmax * dmbx;
              syy += dmay * dmby;
              szz += dmaz * dmbz;
              sxy += dmax * dmby;
              sxz += dmax * dmbz;
              syz += dmay * dmbz;
              ++n;
            }
            const double inv = (n > 0 && std::isfinite(denom_tensor)) ? (1.0 / (static_cast<double>(n) * denom_tensor)) : 0.0;
            const double lxx = sxx * inv;
            const double lyy = syy * inv;
            const double lzz = szz * inv;
            const double lxy = sxy * inv;
            const double lxz = sxz * inv;
            const double lyz = syz * inv;
            ofs << lag << ' '
                << std::setprecision(17) << dt << ' '
                << species_.species_id[a] << ' '
                << species_.species_id[b] << ' '
                << std::setprecision(17) << lxx << ' '
                << std::setprecision(17) << lyy << ' '
                << std::setprecision(17) << lzz << ' '
                << std::setprecision(17) << lxy << ' '
                << std::setprecision(17) << lxz << ' '
                << std::setprecision(17) << lyz << ' '
                << std::setprecision(17) << (lxx + lyy + lzz) / 3.0 << ' '
                << n << ' '
                << std::setprecision(17) << dt << '\n';
          }
        }
      }
    });

    util::atomic_write_text(transference_output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: transference-number-like row sums derived from the Onsager matrix\n";
      ofs << "# definition: t_i(t) = sum_j L_ij^trace / sum_{ab} L_ab^trace (lab or drift-corrected frame, depending on remove_drift)\n";
      ofs << "# columns: lag time species_id t_plus_like row_sum_trace total_trace count_pairs mean_dtimestep\n";
      for (const std::size_t lag : lags) {
        const double dt = mean_dtimestep_for_lag(timesteps_, lag, static_cast<std::size_t>(opt_.lag.origin_stride));
        const double denom_tensor = (dt > 0.0 && mean_volume > 0.0) ? (2.0 * mean_volume * opt_.kbt * dt) : std::numeric_limits<double>::quiet_NaN();
        std::vector<double> row_sum(species_.species_id.size(), 0.0);
        std::size_t n_pairs = 0;
        for (std::size_t a = 0; a < species_.species_id.size(); ++a) {
          for (std::size_t b = 0; b < species_.species_id.size(); ++b) {
            double strace = 0.0;
            std::size_t n = 0;
            for (std::size_t o = 0; o + lag < x_series_.size(); o += static_cast<std::size_t>(opt_.lag.origin_stride)) {
              double dmax = 0.0, dmay = 0.0, dmaz = 0.0;
              double dmbx = 0.0, dmby = 0.0, dmbz = 0.0;
              for (const std::size_t p : species_.selpos_by_species[a]) {
                const double qi = q_sel_[p];
                dmax += qi * (x_series_[o + lag][p] - x_series_[o][p]);
                dmay += qi * (y_series_[o + lag][p] - y_series_[o][p]);
                dmaz += qi * (z_series_[o + lag][p] - z_series_[o][p]);
              }
              for (const std::size_t p : species_.selpos_by_species[b]) {
                const double qi = q_sel_[p];
                dmbx += qi * (x_series_[o + lag][p] - x_series_[o][p]);
                dmby += qi * (y_series_[o + lag][p] - y_series_[o][p]);
                dmbz += qi * (z_series_[o + lag][p] - z_series_[o][p]);
              }
              strace += dmax * dmbx + dmay * dmby + dmaz * dmbz;
              ++n;
            }
            if (n > 0 && std::isfinite(denom_tensor)) {
              row_sum[a] += strace / (static_cast<double>(n) * 3.0 * denom_tensor);
              n_pairs = n;
            }
          }
        }
        const double total = std::accumulate(row_sum.begin(), row_sum.end(), 0.0);
        for (std::size_t a = 0; a < species_.species_id.size(); ++a) {
          const double tlike = (std::abs(total) > 0.0) ? (row_sum[a] / total) : 0.0;
          ofs << lag << ' '
              << std::setprecision(17) << dt << ' '
              << species_.species_id[a] << ' '
              << std::setprecision(17) << tlike << ' '
              << std::setprecision(17) << row_sum[a] << ' '
              << std::setprecision(17) << total << ' '
              << n_pairs << ' '
              << std::setprecision(17) << dt << '\n';
        }
      }
    });
  }

  std::string type_name_;
  std::string instance_name_;
  std::string matrix_output_path_;
  std::string transference_output_path_;
  std::string sel_name_owned_;
  std::string drift_name_owned_;
  SelectionView sel_;
  SelectionView drift_sel_;
  SpeciesGroups species_;
  Options opt_;
  bool started_ = false;
  std::vector<double> q_sel_;
  std::vector<std::vector<double>> x_series_, y_series_, z_series_;
  std::vector<std::int64_t> timesteps_;
  double volume_sum_ = 0.0;
};

void append_species_transport_caps(const IniConfig& cfg,
                                   const std::string& section,
                                   MeasureCapabilities& caps) {
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_identity_consistent = true;
  caps.scale = ScaleCompatibility{true, true, true};
  caps.group_refs.push_back(cfg.get_string(section, "group", std::optional<std::string>("all")));
  if (cfg.has_key(section, "drift_group")) caps.group_refs.push_back(cfg.get_string(section, "drift_group"));
  const std::string species_field = cfg.get_string(section, "species_field", std::optional<std::string>("type"));
  if (species_field == "type") {
    caps.requires_intfields.push_back("type");
  } else if (species_field == "id" || species_field == "mol") {
    caps.requires_i64fields.push_back(species_field);
  } else {
    caps.requires_dfields.push_back(species_field);
  }
}

std::unique_ptr<IMeasure> create_species_tensor(const IniConfig& cfg,
                                                const std::string& section,
                                                const std::string& instance,
                                                const MeasureBuildEnv& env,
                                                const SystemContext& sysctx,
                                                SpeciesTensorMode mode,
                                                std::string type_name) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error(type_name + " factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, type_name);
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(false));
  const std::string drift_group = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));
  SelectionView drift_sel = remove_drift ? get_static_group_view(*env.selection_provider, frame0, drift_group, type_name)
                                         : get_static_group_view(*env.selection_provider, frame0, "all", type_name);

  const std::string species_field = cfg.get_string(section, "species_field", std::optional<std::string>("type"));
  const auto species_per_atom = integer_like_field_to_i64(frame0, species_field, true);
  SpeciesGroups species = build_species_groups(sel, species_per_atom);

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const std::string default_out = (mode == SpeciesTensorMode::MSD) ? (type_name + std::string(".dat")) : "vacf.dat";
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>(default_out))).lexically_normal();

  SpeciesTensorOptions opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = resolve_exact_frame_end(parse_correlator_spec(cfg, section, env.dt), env.follow,
                                                opt.range.frame_start,
                                                cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1)),
                                                env.input_path);
  opt.range.dry_run = env.dry_run;
  opt.corr = parse_correlator_spec(cfg, section, env.dt);
  opt.mode = mode;
  opt.remove_drift = remove_drift;
  opt.vector_x_field = cfg.get_string(section, "vector_x_field", std::optional<std::string>("vx"));
  opt.vector_y_field = cfg.get_string(section, "vector_y_field", std::optional<std::string>("vy"));
  opt.vector_z_field = cfg.get_string(section, "vector_z_field", std::optional<std::string>("vz"));
  return std::make_unique<SpeciesTensorCorrMeasure>(type_name, instance, out.string(), sel, drift_sel, opt, std::move(species));
}

MeasureCapabilities species_msd_caps(const IniConfig& cfg,
                                     const std::string& section,
                                     const std::string& instance,
                                     const MeasureBuildEnv& env) {
  (void)instance; (void)env;
  MeasureCapabilities caps;
  append_species_transport_caps(cfg, section, caps);
  caps.requires_dfields = {"xu", "yu", "zu"};
  return caps;
}

std::unique_ptr<IMeasure> species_msd_create(const IniConfig& cfg,
                                             const std::string& section,
                                             const std::string& instance,
                                             const MeasureBuildEnv& env,
                                             const SystemContext& sysctx) {
  return create_species_tensor(cfg, section, instance, env, sysctx, SpeciesTensorMode::MSD, "species_msd");
}

std::unique_ptr<IMeasure> diffusion_tensor_create(const IniConfig& cfg,
                                                  const std::string& section,
                                                  const std::string& instance,
                                                  const MeasureBuildEnv& env,
                                                  const SystemContext& sysctx) {
  return create_species_tensor(cfg, section, instance, env, sysctx, SpeciesTensorMode::MSD, "diffusion_tensor");
}

MeasureCapabilities vacf_caps(const IniConfig& cfg,
                              const std::string& section,
                              const std::string& instance,
                              const MeasureBuildEnv& env) {
  (void)instance; (void)env;
  MeasureCapabilities caps;
  append_species_transport_caps(cfg, section, caps);
  caps.requires_dfields = {cfg.get_string(section, "vector_x_field", std::optional<std::string>("vx")),
                           cfg.get_string(section, "vector_y_field", std::optional<std::string>("vy")),
                           cfg.get_string(section, "vector_z_field", std::optional<std::string>("vz"))};
  return caps;
}

std::unique_ptr<IMeasure> vacf_create(const IniConfig& cfg,
                                      const std::string& section,
                                      const std::string& instance,
                                      const MeasureBuildEnv& env,
                                      const SystemContext& sysctx) {
  return create_species_tensor(cfg, section, instance, env, sysctx, SpeciesTensorMode::VACF, "vacf");
}

MeasureCapabilities species_current_corr_caps(const IniConfig& cfg,
                                              const std::string& section,
                                              const std::string& instance,
                                              const MeasureBuildEnv& env) {
  (void)instance; (void)env;
  MeasureCapabilities caps;
  append_species_transport_caps(cfg, section, caps);
  caps.requires_dfields = {cfg.get_string(section, "q_field", std::optional<std::string>("q")),
                           cfg.get_string(section, "vector_x_field", std::optional<std::string>("vx")),
                           cfg.get_string(section, "vector_y_field", std::optional<std::string>("vy")),
                           cfg.get_string(section, "vector_z_field", std::optional<std::string>("vz"))};
  return caps;
}

std::unique_ptr<IMeasure> species_current_corr_create(const IniConfig& cfg,
                                                      const std::string& section,
                                                      const std::string& instance,
                                                      const MeasureBuildEnv& env,
                                                      const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("species_current_corr factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "species_current_corr");
  const std::string species_field = cfg.get_string(section, "species_field", std::optional<std::string>("type"));
  const auto species_per_atom = integer_like_field_to_i64(frame0, species_field, true);
  SpeciesGroups species = build_species_groups(sel, species_per_atom);

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("species_current_corr.dat"))).lexically_normal();

  SpeciesCurrentCorrMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.lag = lag_options_from_config(cfg, section);
  opt.q_field = cfg.get_string(section, "q_field", std::optional<std::string>("q"));
  opt.vx_field = cfg.get_string(section, "vector_x_field", std::optional<std::string>("vx"));
  opt.vy_field = cfg.get_string(section, "vector_y_field", std::optional<std::string>("vy"));
  opt.vz_field = cfg.get_string(section, "vector_z_field", std::optional<std::string>("vz"));
  return std::make_unique<SpeciesCurrentCorrMeasure>(instance, out.string(), sel, std::move(species), opt);
}

MeasureCapabilities onsager_transport_caps(const IniConfig& cfg,
                                           const std::string& section,
                                           const std::string& instance,
                                           const MeasureBuildEnv& env) {
  (void)instance; (void)env;
  MeasureCapabilities caps;
  append_species_transport_caps(cfg, section, caps);
  caps.requires_dfields = {"xu", "yu", "zu", cfg.get_string(section, "q_field", std::optional<std::string>("q"))};
  return caps;
}

std::unique_ptr<IMeasure> create_onsager_like(const IniConfig& cfg,
                                              const std::string& section,
                                              const std::string& instance,
                                              const MeasureBuildEnv& env,
                                              const SystemContext& sysctx,
                                              std::string type_name) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error(type_name + " factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, type_name);
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(false));
  const std::string drift_group = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));
  SelectionView drift_sel = remove_drift ? get_static_group_view(*env.selection_provider, frame0, drift_group, type_name)
                                         : get_static_group_view(*env.selection_provider, frame0, "all", type_name);
  const std::string species_field = cfg.get_string(section, "species_field", std::optional<std::string>("type"));
  const auto species_per_atom = integer_like_field_to_i64(frame0, species_field, true);
  SpeciesGroups species = build_species_groups(sel, species_per_atom);

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out_matrix = (output_dir / cfg.get_string(section, "output_matrix", std::optional<std::string>(type_name + std::string("_matrix.dat")))).lexically_normal();
  const fs::path out_t = (output_dir / cfg.get_string(section, "output_transference", std::optional<std::string>(type_name + std::string("_transference.dat")))).lexically_normal();

  OnsagerTransportMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.lag = lag_options_from_config(cfg, section);
  opt.q_field = cfg.get_string(section, "q_field", std::optional<std::string>("q"));
  opt.remove_drift = remove_drift;
  opt.kbt = cfg.get_double(section, "kbt");
  return std::make_unique<OnsagerTransportMeasure>(type_name, instance, out_matrix.string(), out_t.string(), sel, drift_sel, std::move(species), opt);
}

std::unique_ptr<IMeasure> onsager_transport_create(const IniConfig& cfg,
                                                   const std::string& section,
                                                   const std::string& instance,
                                                   const MeasureBuildEnv& env,
                                                   const SystemContext& sysctx) {
  return create_onsager_like(cfg, section, instance, env, sysctx, "onsager_transport");
}

std::unique_ptr<IMeasure> transference_number_create(const IniConfig& cfg,
                                                     const std::string& section,
                                                     const std::string& instance,
                                                     const MeasureBuildEnv& env,
                                                     const SystemContext& sysctx) {
  return create_onsager_like(cfg, section, instance, env, sysctx, "transference_number");
}

static MeasureRegistrar g_reg_species_msd("species_msd", &species_msd_caps, &species_msd_create);
static MeasureRegistrar g_reg_diffusion_tensor("diffusion_tensor", &species_msd_caps, &diffusion_tensor_create);
static MeasureRegistrar g_reg_vacf("vacf", &vacf_caps, &vacf_create);
static MeasureRegistrar g_reg_species_current_corr("species_current_corr", &species_current_corr_caps, &species_current_corr_create);
static MeasureRegistrar g_reg_onsager_transport("onsager_transport", &onsager_transport_caps, &onsager_transport_create);
static MeasureRegistrar g_reg_transference_number("transference_number", &onsager_transport_caps, &transference_number_create);

} // namespace
} // namespace pilots
