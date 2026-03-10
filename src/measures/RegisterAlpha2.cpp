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
using measure_ext::parse_diag_mask;
using measure_ext::resolve_exact_frame_end;
using measure_ext::resolve_path;
using measure_ext::x_unit_for_axis;

struct Alpha2PairOpT6 {
  int diag_mask = 7;

  Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const {
    const bool use_x = (diag_mask & 1) != 0;
    const bool use_y = (diag_mask & 2) != 0;
    const bool use_z = (diag_mask & 4) != 0;

    const std::size_t n = cur.xx.size();
    if (org.xx.size() != n || cur.yy.size() != n || org.yy.size() != n ||
        cur.zz.size() != n || org.zz.size() != n) {
      throw std::runtime_error("Alpha2PairOpT6: slot sizes mismatch");
    }

    double sum_r2 = 0.0;
    double sum_r4 = 0.0;
#if PILOTS_HAS_OPENMP
#pragma omp parallel for reduction(+:sum_r2,sum_r4)
#endif
    for (std::size_t p = 0; p < n; ++p) {
      const double dx = cur.xx[p] - org.xx[p];
      const double dy = cur.yy[p] - org.yy[p];
      const double dz = cur.zz[p] - org.zz[p];
      double r2 = 0.0;
      if (use_x) r2 += dx * dx;
      if (use_y) r2 += dy * dy;
      if (use_z) r2 += dz * dz;
      sum_r2 += r2;
      sum_r4 += r2 * r2;
    }

    Tensor6 out;
    out.v[T6_XX] = sum_r2 / static_cast<double>(n);
    out.v[T6_YY] = sum_r4 / static_cast<double>(n);
    return out;
  }
};

class Alpha2Measure final : public IMeasure {
public:
  struct Options {
    std::int64_t frame_start = 0;
    std::int64_t frame_end = -1;
    int diag_mask = 7;
    bool remove_drift = true;
    bool dry_run = false;
    CorrelatorSpec corr;
  };

  Alpha2Measure(std::string instance_name,
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

    if (opt_.frame_start < 0) throw std::runtime_error("Alpha2Measure: frame_start must be >= 0");
    if (opt_.frame_end >= 0 && opt_.frame_end < opt_.frame_start) {
      throw std::runtime_error("Alpha2Measure: frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.corr.dt > 0.0)) throw std::runtime_error("Alpha2Measure: corr.dt must be > 0");
    if (opt_.diag_mask <= 0 || opt_.diag_mask > 7) throw std::runtime_error("Alpha2Measure: invalid diag_mask");
    if (sel_.idx.empty()) throw std::runtime_error("Alpha2Measure: selection is empty");
    if (opt_.remove_drift && drift_sel_.idx.empty()) {
      throw std::runtime_error("Alpha2Measure: drift selection is empty");
    }
    if (opt_.corr.axis == LagAxis::TimeBin && !(opt_.corr.timebin_width > 0.0)) {
      throw std::runtime_error("Alpha2Measure: timebin_width must be > 0 for lag_axis=timebin");
    }
    if (opt_.corr.type == "exact" && opt_.frame_end < 0) {
      throw std::runtime_error("Alpha2Measure: correlator=exact requires finite frame_end");
    }

    std::size_t window_frames = 0;
    if (opt_.corr.type == "exact") {
      window_frames = static_cast<std::size_t>(opt_.frame_end - opt_.frame_start + 1);
      if (window_frames < 2) {
        throw std::runtime_error("Alpha2Measure: exact window must have >= 2 frames");
      }
    }

    const Alpha2PairOpT6 pair_op{opt_.diag_mask};
    corr_ = make_correlator_t6_runtime_auto(sel_.idx.size(), window_frames, opt_.corr, pair_op);
  }

  std::string type() const override { return "alpha2"; }
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
    od.columns = {"lag", "time", "r2", "r4", "alpha2", "count_pairs", "sem_r2", "sem_r4", "n_blocks", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));

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
    if (!opt_.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    const std::int64_t fi = static_cast<std::int64_t>(frame_index);
    if (fi < opt_.frame_start) return;
    if (opt_.frame_end >= 0 && fi > opt_.frame_end) return;

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

    corr_->push(frame.timestep,
                std::span<const double>(tmp_xx_.data(), tmp_xx_.size()),
                std::span<const double>(tmp_yy_.data(), tmp_yy_.size()),
                std::span<const double>(tmp_zz_.data(), tmp_zz_.size()),
                std::span<const double>(tmp_xy_.data(), tmp_xy_.size()),
                std::span<const double>(tmp_xz_.data(), tmp_xz_.size()),
                std::span<const double>(tmp_yz_.data(), tmp_yz_.size()));
  }

  void flush_partial() override {
    if (opt_.dry_run || !started_) return;
    const CorrelationSeriesT6 series = corr_->snapshot();
    const int dim = count_diag_dims(opt_.diag_mask);
    const double pref = static_cast<double>(dim) / static_cast<double>(dim + 2);

    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      write_header_(ofs);
      const std::size_t m = series.lag.size();
      for (std::size_t i = 0; i < m; ++i) {
        const double r2 = series.value_xx[i];
        const double r4 = series.value_yy[i];
        const double denom = r2 * r2;
        const double a2 = (denom > 0.0) ? (pref * r4 / denom - 1.0) : 0.0;
        ofs << std::setprecision(17) << series.lag[i] << ' '
            << std::setprecision(17) << series.time[i] << ' '
            << std::setprecision(17) << r2 << ' '
            << std::setprecision(17) << r4 << ' '
            << std::setprecision(17) << a2 << ' '
            << series.count_pairs[i] << ' '
            << std::setprecision(17) << series.sem_xx[i] << ' '
            << std::setprecision(17) << series.sem_yy[i] << ' '
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
  Options opt_;
  std::unique_ptr<ICorrelatorT6> corr_;
  bool started_ = false;

  std::vector<double> tmp_xx_, tmp_yy_, tmp_zz_, tmp_xy_, tmp_xz_, tmp_yz_;

  void write_header_(std::ostream& ofs) const {
    ofs << "# PILOTS: alpha2(t) using T6 correlator channels XX=<r2>, YY=<r4>\n";
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
    if (opt_.corr.type == "exact") {
      ofs << "# exact_lag_stride: " << opt_.corr.lag_stride << "\n";
    } else if (opt_.corr.type == "multitau") {
      ofs << "# multitau_channels: " << opt_.corr.mt_channels << ", levels: " << opt_.corr.mt_levels << "\n";
    }
    ofs << "# block_size: " << opt_.corr.block_size << " (0 disables SEM)\n";
    ofs << "# columns: lag  time  r2  r4  alpha2  count_pairs  sem_r2  sem_r4  n_blocks  mean_dtimestep\n";
  }
};

MeasureCapabilities alpha2_caps(const IniConfig& cfg,
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

std::unique_ptr<IMeasure> alpha2_create(const IniConfig& cfg,
                                        const std::string& section,
                                        const std::string& instance,
                                        const MeasureBuildEnv& env,
                                        const SystemContext& sysctx) {
  (void)sysctx;

  if (!env.first_frame) throw std::runtime_error("alpha2 factory: first_frame is null");
  if (!env.selection_provider) {
    throw std::runtime_error("alpha2 factory: SelectionProvider is missing");
  }
  const Frame& frame0 = *env.first_frame;

  const std::string group_ref = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo_group_ref = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string combine_expr = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  const bool remove_drift = cfg.get_bool(section, "remove_drift", std::optional<bool>(true));
  const std::string drift_group_ref = cfg.get_string(section, "drift_group", std::optional<std::string>("all"));

  const CorrelatorSpec corr_spec = parse_correlator_spec(cfg, section, env.dt);
  const std::int64_t frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  const std::int64_t frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  const std::string comp_s = cfg.get_string(section, "components", std::optional<std::string>("xxyyzz"));

  if (frame_start < 0) throw std::runtime_error("frame_start must be >= 0");
  if (frame_end >= 0 && frame_end < frame_start) {
    throw std::runtime_error("frame_end must be -1 or >= frame_start");
  }

  const fs::path out_path = measure_ext::resolve_measure_output_path(cfg, section, env, "output", "alpha2.dat");

  const std::int64_t frame_end_eff = resolve_exact_frame_end(corr_spec, env.follow, frame_start, frame_end, env.input_path);
  const int diag_mask = parse_diag_mask(comp_s);

  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group_ref, topo_group_ref, combine_expr, "alpha2");
  SelectionView drift_sel = remove_drift
      ? get_static_group_view(*env.selection_provider, frame0, drift_group_ref, "alpha2 drift_group")
      : get_static_group_view(*env.selection_provider, frame0, "all", "alpha2 drift_group");

  Alpha2Measure::Options opt;
  opt.frame_start = frame_start;
  opt.frame_end = (corr_spec.type == "exact") ? frame_end_eff : frame_end;
  opt.diag_mask = diag_mask;
  opt.remove_drift = remove_drift;
  opt.dry_run = env.dry_run;
  opt.corr = corr_spec;

  return std::make_unique<Alpha2Measure>(instance, out_path.string(), sel, drift_sel, opt);
}

static MeasureRegistrar g_register_alpha2("alpha2", &alpha2_caps, &alpha2_create);

} // namespace
} // namespace pilots
