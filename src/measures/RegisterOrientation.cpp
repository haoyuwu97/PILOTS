#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
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

using measure_ext::OrderedChains;
using measure_ext::SelectedChains;
using measure_ext::Vec3;
using measure_ext::atom_vec3;
using measure_ext::build_ordered_chains_from_chain_pos;
using measure_ext::build_selected_chains;
using measure_ext::chain_id_per_atom_from_config;
using measure_ext::cross;
using measure_ext::dot;
using measure_ext::entity_id_per_atom_from_config;
using measure_ext::get_static_combined_view;
using measure_ext::norm;
using measure_ext::ordered_chains_from_config;
using measure_ext::resolve_exact_frame_end;
using measure_ext::resolve_path;
using measure_ext::safe_unit;
using measure_ext::x_unit_for_axis;

inline double legendre_p(int l, double x) {
  x = std::clamp(x, -1.0, 1.0);
  switch (l) {
    case 0: return 1.0;
    case 1: return x;
    case 2: return 0.5 * (3.0 * x * x - 1.0);
    case 3: return 0.5 * (5.0 * x * x * x - 3.0 * x);
    default:
      throw std::runtime_error("orientation ACF currently supports legendre_order = 0,1,2,3");
  }
}

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

struct LegendrePairOpT6 {
  int legendre_order = 2;

  Tensor6 operator()(const T6Slot& cur, const T6Slot& org) const {
    const std::size_t n = cur.xx.size();
    if (org.xx.size() != n || cur.yy.size() != n || org.yy.size() != n ||
        cur.zz.size() != n || org.zz.size() != n) {
      throw std::runtime_error("LegendrePairOpT6: slot sizes mismatch");
    }
    double sum = 0.0;
#if PILOTS_HAS_OPENMP
#pragma omp parallel for reduction(+:sum)
#endif
    for (std::size_t p = 0; p < n; ++p) {
      const double d = cur.xx[p] * org.xx[p] + cur.yy[p] * org.yy[p] + cur.zz[p] * org.zz[p];
      sum += legendre_p(legendre_order, d);
    }
    Tensor6 out;
    out.v[T6_XX] = sum / static_cast<double>(n);
    return out;
  }
};

class OrientationCorrMeasure final : public IMeasure {
public:
  using VectorBuilderFn = std::function<void(const Frame&, std::vector<double>&, std::vector<double>&, std::vector<double>&)>;

  struct Options {
    RangeOptions range;
    CorrelatorSpec corr;
    int legendre_order = 2;
    std::string object_kind;
  };

  OrientationCorrMeasure(std::string type_name,
                         std::string instance_name,
                         std::string output_path,
                         SelectionView sel,
                         std::size_t n_objects,
                         Options opt,
                         VectorBuilderFn builder)
      : type_name_(std::move(type_name)),
        instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        opt_(std::move(opt)),
        builder_(std::move(builder)) {
    sel_name_owned_ = std::string(sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    if (sel_.idx.empty()) throw std::runtime_error(type_name_ + ": selection is empty");
    if (n_objects == 0) throw std::runtime_error(type_name_ + ": no objects were constructed");
    if (opt_.range.frame_start < 0) throw std::runtime_error(type_name_ + ": frame_start must be >= 0");
    if (opt_.range.frame_end >= 0 && opt_.range.frame_end < opt_.range.frame_start) {
      throw std::runtime_error(type_name_ + ": frame_end must be -1 or >= frame_start");
    }
    corr_ = make_correlator_t6_runtime_auto(n_objects, exact_window_frames_(), opt_.corr,
                                            LegendrePairOpT6{opt_.legendre_order});
    tmp_xx_.assign(n_objects, 0.0);
    tmp_yy_.assign(n_objects, 0.0);
    tmp_zz_.assign(n_objects, 0.0);
    tmp_zero_.assign(n_objects, 0.0);
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
    od.columns = {"lag", "time", "c_l", "c_l_norm", "count_pairs", "sem", "n_blocks", "mean_dtimestep"};
    md.outputs.push_back(std::move(od));

    md.params["legendre_order"] = std::to_string(opt_.legendre_order);
    md.params["object_kind"] = opt_.object_kind;
    md.params["frame_start"] = std::to_string(opt_.range.frame_start);
    md.params["frame_end"] = std::to_string(opt_.range.frame_end);
    md.params["correlator"] = opt_.corr.type;
    md.params["lag_axis"] = lag_axis_name(opt_.corr.axis);
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
    builder_(frame, tmp_xx_, tmp_yy_, tmp_zz_);
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
    const CorrelationSeriesT6 series = corr_->snapshot();
    const double c0 = !series.value_xx.empty() ? series.value_xx.front() : 0.0;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: orientation / reorientation autocorrelation\n";
      ofs << "# type: " << type_name_ << ", object_kind: " << opt_.object_kind << "\n";
      ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
      ofs << "# legendre_order: " << opt_.legendre_order << "\n";
      ofs << "# frame_range: [" << opt_.range.frame_start << ", " << opt_.range.frame_end << "]\n";
      ofs << "# correlator: " << opt_.corr.type << ", lag_axis: " << lag_axis_name(opt_.corr.axis) << "\n";
      ofs << "# columns: lag time c_l c_l_norm count_pairs sem n_blocks mean_dtimestep\n";
      for (std::size_t i = 0; i < series.lag.size(); ++i) {
        const double c = series.value_xx[i];
        const double cn = (std::abs(c0) > 0.0) ? (c / c0) : 0.0;
        ofs << std::setprecision(17) << series.lag[i] << ' '
            << std::setprecision(17) << series.time[i] << ' '
            << std::setprecision(17) << c << ' '
            << std::setprecision(17) << cn << ' '
            << series.count_pairs[i] << ' '
            << std::setprecision(17) << series.sem_xx[i] << ' '
            << series.n_blocks[i] << ' '
            << std::setprecision(17) << series.mean_dtimestep[i] << '\n';
      }
    });
  }

  void finalize() override {
    if (opt_.range.dry_run) return;
    flush_partial();
  }

private:
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
  SelectionView sel_;
  Options opt_;
  VectorBuilderFn builder_;
  std::unique_ptr<ICorrelatorT6> corr_;
  bool started_ = false;
  std::vector<double> tmp_xx_, tmp_yy_, tmp_zz_, tmp_zero_;
};

class TotalDipoleTimeMeasure final : public IMeasure {
public:
  struct Options {
    RangeOptions range;
    std::int64_t max_lag_frames = -1;
    std::int64_t lag_stride = 1;
    std::int64_t origin_stride = 1;
    std::string charge_field = "q";
    double delta_epsilon = std::numeric_limits<double>::quiet_NaN();
    double epsilon_infty = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> omega_list;
    bool write_spectrum = false;
  };

  TotalDipoleTimeMeasure(std::string type_name,
                         std::string instance_name,
                         std::string acf_output_path,
                         std::string spectrum_output_path,
                         SelectionView sel,
                         SelectedChains entities,
                         Options opt)
      : type_name_(std::move(type_name)),
        instance_name_(std::move(instance_name)),
        acf_output_path_(std::move(acf_output_path)),
        spectrum_output_path_(std::move(spectrum_output_path)),
        entities_(std::move(entities)),
        opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    if (sel_.idx.empty()) throw std::runtime_error(type_name_ + ": selection is empty");
    if (entities_.n_chains() == 0) throw std::runtime_error(type_name_ + ": no entities were constructed");
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
    od1.path = acf_output_path_;
    od1.format = "text";
    od1.x_axis = "lag";
    od1.x_unit = "frames";
    od1.columns = {"lag", "time", "c_xx", "c_yy", "c_zz", "c_xy", "c_xz", "c_yz", "c_trace", "c_trace_norm", "count_pairs", "mean_dtimestep"};
    md.outputs.push_back(std::move(od1));

    if (opt_.write_spectrum) {
      output::OutputFileDescriptor od2;
      od2.path = spectrum_output_path_;
      od2.format = "text";
      od2.x_axis = "omega";
      od2.x_unit = "inverse_time";
      od2.columns = {"omega", "phi_storage", "phi_loss", "epsilon_prime", "epsilon_double_prime"};
      md.outputs.push_back(std::move(od2));
    }

    md.params["charge_field"] = opt_.charge_field;
    return md;
  }

  void on_start(const Frame& first_frame) override {
    const auto q = first_frame.require_dfield(opt_.charge_field);
    q_sel_.assign(sel_.idx.size(), 0.0);
    for (std::size_t p = 0; p < sel_.idx.size(); ++p) q_sel_[p] = q[sel_.idx[p]];
    started_ = true;
    if (!opt_.range.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.range)) return;
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");
    Vec3 mtot{};
    for (std::size_t c = 0; c < entities_.n_chains(); ++c) {
      Vec3 com{};
      const double inv = 1.0 / static_cast<double>(entities_.selpos_by_chain[c].size());
      for (const std::size_t p : entities_.selpos_by_chain[c]) {
        const std::size_t atom = sel_.idx[p];
        com.x += xu[atom];
        com.y += yu[atom];
        com.z += zu[atom];
      }
      com = inv * com;
      Vec3 mu{};
      for (const std::size_t p : entities_.selpos_by_chain[c]) {
        const std::size_t atom = sel_.idx[p];
        const double qi = q_sel_[p];
        mu.x += qi * (xu[atom] - com.x);
        mu.y += qi * (yu[atom] - com.y);
        mu.z += qi * (zu[atom] - com.z);
      }
      mtot += mu;
    }
    m_series_.push_back(mtot);
    timesteps_.push_back(frame.timestep);
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
    if (m_series_.empty()) return;
    const std::size_t nframes = m_series_.size();
    const std::size_t max_lag = (opt_.max_lag_frames < 0)
        ? (nframes - 1)
        : std::min<std::size_t>(static_cast<std::size_t>(opt_.max_lag_frames), nframes - 1);
    std::vector<std::size_t> lags;
    for (std::size_t lag = 0; lag <= max_lag; lag += static_cast<std::size_t>(opt_.lag_stride)) lags.push_back(lag);
    if (lags.empty()) lags.push_back(0);

    struct Row {
      std::size_t lag = 0;
      double time = 0.0;
      double cxx = 0.0, cyy = 0.0, czz = 0.0, cxy = 0.0, cxz = 0.0, cyz = 0.0;
      std::size_t count = 0;
    };
    std::vector<Row> rows;
    rows.reserve(lags.size());
    for (const std::size_t lag : lags) {
      Row r;
      r.lag = lag;
      r.time = 0.0;
      for (std::size_t o = 0; o + lag < nframes; o += static_cast<std::size_t>(opt_.origin_stride)) {
        const Vec3& a = m_series_[o + lag];
        const Vec3& b = m_series_[o];
        r.cxx += a.x * b.x;
        r.cyy += a.y * b.y;
        r.czz += a.z * b.z;
        r.cxy += a.x * b.y;
        r.cxz += a.x * b.z;
        r.cyz += a.y * b.z;
        r.time += static_cast<double>(timesteps_[o + lag] - timesteps_[o]);
        ++r.count;
      }
      if (r.count > 0) {
        const double inv = 1.0 / static_cast<double>(r.count);
        r.cxx *= inv; r.cyy *= inv; r.czz *= inv; r.cxy *= inv; r.cxz *= inv; r.cyz *= inv; r.time *= inv;
      }
      rows.push_back(r);
    }

    const double c0 = rows.empty() ? 0.0 : (rows.front().cxx + rows.front().cyy + rows.front().czz);
    util::atomic_write_text(acf_output_path_, [&](std::ostream& ofs) {
      ofs << "# PILOTS: total-dipole autocorrelation\n";
      ofs << "# type: " << type_name_ << "\n";
      ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
      ofs << "# charge_field: " << opt_.charge_field << "\n";
      ofs << "# columns: lag time c_xx c_yy c_zz c_xy c_xz c_yz c_trace c_trace_norm count_pairs mean_dtimestep\n";
      for (const auto& r : rows) {
        const double trace = r.cxx + r.cyy + r.czz;
        const double tn = (std::abs(c0) > 0.0) ? (trace / c0) : 0.0;
        ofs << r.lag << ' '
            << std::setprecision(17) << r.time << ' '
            << std::setprecision(17) << r.cxx << ' '
            << std::setprecision(17) << r.cyy << ' '
            << std::setprecision(17) << r.czz << ' '
            << std::setprecision(17) << r.cxy << ' '
            << std::setprecision(17) << r.cxz << ' '
            << std::setprecision(17) << r.cyz << ' '
            << std::setprecision(17) << trace << ' '
            << std::setprecision(17) << tn << ' '
            << r.count << ' '
            << std::setprecision(17) << r.time << '\n';
      }
    });

    if (opt_.write_spectrum && !opt_.omega_list.empty()) {
      util::atomic_write_text(spectrum_output_path_, [&](std::ostream& ofs) {
        ofs << "# PILOTS: dielectric spectrum from normalized dipole relaxation function\n";
        ofs << "# formula: phi_storage = 1 - omega * integral Phi(t) sin(omega t) dt; phi_loss = omega * integral Phi(t) cos(omega t) dt\n";
        ofs << "# columns: omega phi_storage phi_loss epsilon_prime epsilon_double_prime\n";
        const std::vector<double> t(rows.size(), 0.0);
        std::vector<double> phi(rows.size(), 0.0);
        for (std::size_t i = 0; i < rows.size(); ++i) {
          const double trace = rows[i].cxx + rows[i].cyy + rows[i].czz;
          phi[i] = (std::abs(c0) > 0.0) ? (trace / c0) : 0.0;
        }
        std::vector<double> tt(rows.size(), 0.0);
        for (std::size_t i = 0; i < rows.size(); ++i) tt[i] = rows[i].time;
        for (const double omega : opt_.omega_list) {
          double int_s = 0.0, int_c = 0.0;
          for (std::size_t i = 1; i < rows.size(); ++i) {
            const double dt = tt[i] - tt[i - 1];
            const double y1s = phi[i - 1] * std::sin(omega * tt[i - 1]);
            const double y2s = phi[i] * std::sin(omega * tt[i]);
            const double y1c = phi[i - 1] * std::cos(omega * tt[i - 1]);
            const double y2c = phi[i] * std::cos(omega * tt[i]);
            int_s += 0.5 * dt * (y1s + y2s);
            int_c += 0.5 * dt * (y1c + y2c);
          }
          const double phi_storage = 1.0 - omega * int_s;
          const double phi_loss = omega * int_c;
          double epr = std::numeric_limits<double>::quiet_NaN();
          double epi = std::numeric_limits<double>::quiet_NaN();
          if (std::isfinite(opt_.delta_epsilon) && std::isfinite(opt_.epsilon_infty)) {
            epr = opt_.epsilon_infty + opt_.delta_epsilon * phi_storage;
            epi = opt_.delta_epsilon * phi_loss;
          }
          ofs << std::setprecision(17) << omega << ' '
              << std::setprecision(17) << phi_storage << ' '
              << std::setprecision(17) << phi_loss << ' '
              << std::setprecision(17) << epr << ' '
              << std::setprecision(17) << epi << '\n';
        }
      });
    }
  }

  std::string type_name_;
  std::string instance_name_;
  std::string acf_output_path_;
  std::string spectrum_output_path_;
  std::string sel_name_owned_;
  SelectionView sel_;
  SelectedChains entities_;
  Options opt_;
  bool started_ = false;
  std::vector<double> q_sel_;
  std::vector<Vec3> m_series_;
  std::vector<std::int64_t> timesteps_;
};

MeasureCapabilities bond_vector_acf_caps(const IniConfig& cfg,
                                         const std::string& section,
                                         const std::string& instance,
                                         const MeasureBuildEnv& env) {
  (void)instance; (void)env;
  MeasureCapabilities caps;
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_identity_consistent = true;
  caps.requires_dfields = {"xu", "yu", "zu"};
  caps.requires_topology_sections = {"bonds"};
  caps.scale = ScaleCompatibility{true, true, true};
  caps.group_refs.push_back(cfg.get_string(section, "group", std::optional<std::string>("all")));
  return caps;
}

std::unique_ptr<IMeasure> bond_vector_acf_create(const IniConfig& cfg,
                                                 const std::string& section,
                                                 const std::string& instance,
                                                 const MeasureBuildEnv& env,
                                                 const SystemContext& sysctx) {
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("bond_vector_acf factory: missing first_frame or SelectionProvider");
  if (!sysctx.topology || !sysctx.topology->has_bonds()) throw std::runtime_error("bond_vector_acf requires topology bonds");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "bond_vector_acf");

  std::unordered_set<std::size_t> atomset(sel.idx.begin(), sel.idx.end());
  std::vector<int> bond_types;
  if (cfg.has_key(section, "bond_types")) bond_types = measure_ext::parse_int_list(cfg, section, "bond_types");
  std::unordered_set<int> type_set(bond_types.begin(), bond_types.end());
  std::vector<std::pair<std::size_t, std::size_t>> bonds;
  for (const auto& b : sysctx.topology->bonds) {
    if (atomset.count(b.i) == 0 || atomset.count(b.j) == 0) continue;
    if (!type_set.empty() && type_set.count(b.type) == 0) continue;
    bonds.push_back({b.i, b.j});
  }
  if (bonds.empty()) throw std::runtime_error("bond_vector_acf: no bonds remained after selection/filtering");

  const auto builder = [bonds](const Frame& frame, std::vector<double>& xx, std::vector<double>& yy, std::vector<double>& zz) {
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");
#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
    for (std::size_t k = 0; k < bonds.size(); ++k) {
      const auto [i, j] = bonds[k];
      const auto dr = frame.box.min_image_displacement(xu[i], yu[i], zu[i], xu[j], yu[j], zu[j]);
      const Vec3 u = safe_unit(Vec3{dr[0], dr[1], dr[2]});
      xx[k] = u.x; yy[k] = u.y; zz[k] = u.z;
    }
  };

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("bond_vector_acf.dat"))).lexically_normal();

  OrientationCorrMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = resolve_exact_frame_end(parse_correlator_spec(cfg, section, env.dt), env.follow,
                                                opt.range.frame_start,
                                                cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1)),
                                                env.input_path);
  opt.range.dry_run = env.dry_run;
  opt.corr = parse_correlator_spec(cfg, section, env.dt);
  opt.legendre_order = static_cast<int>(cfg.get_int64(section, "legendre_order", std::optional<std::int64_t>(2)));
  opt.object_kind = "bond_vector";
  return std::make_unique<OrientationCorrMeasure>("bond_vector_acf", instance, out.string(), sel, bonds.size(), opt, builder);
}

MeasureCapabilities segment_reorientation_caps(const IniConfig& cfg,
                                               const std::string& section,
                                               const std::string& instance,
                                               const MeasureBuildEnv& env) {
  (void)instance; (void)env;
  MeasureCapabilities caps;
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_identity_consistent = true;
  caps.requires_dfields = {"xu", "yu", "zu"};
  caps.scale = ScaleCompatibility{true, true, true};
  caps.group_refs.push_back(cfg.get_string(section, "group", std::optional<std::string>("all")));
  if (cfg.has_key(section, "chain_id_field")) caps.requires_dfields.push_back(cfg.get_string(section, "chain_id_field"));
  if (cfg.has_key(section, "chain_pos_field")) caps.requires_dfields.push_back(cfg.get_string(section, "chain_pos_field"));
  else caps.requires_topology_sections.push_back("bonds");
  return caps;
}

std::unique_ptr<IMeasure> segment_reorientation_create(const IniConfig& cfg,
                                                       const std::string& section,
                                                       const std::string& instance,
                                                       const MeasureBuildEnv& env,
                                                       const SystemContext& sysctx) {
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("segment_reorientation factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "segment_reorientation");
  const auto chain_ids = chain_id_per_atom_from_config(cfg, section, frame0, sysctx);
  const SelectedChains groups = build_selected_chains(sel, chain_ids);
  const OrderedChains ordered = ordered_chains_from_config(cfg, section, frame0, sysctx, sel, groups);
  const std::size_t span = static_cast<std::size_t>(cfg.get_int64(section, "segment_span", std::optional<std::int64_t>(1)));
  const std::size_t stride = static_cast<std::size_t>(cfg.get_int64(section, "segment_stride", std::optional<std::int64_t>(1)));
  if (span == 0) throw std::runtime_error("segment_reorientation: segment_span must be >= 1");
  if (stride == 0) throw std::runtime_error("segment_reorientation: segment_stride must be >= 1");
  std::vector<std::pair<std::size_t, std::size_t>> segments;
  for (const auto& chain : ordered.atom_by_chain) {
    if (chain.size() <= span) continue;
    for (std::size_t n = 0; n + span < chain.size(); n += stride) {
      segments.push_back({chain[n], chain[n + span]});
    }
  }
  if (segments.empty()) throw std::runtime_error("segment_reorientation: no segments generated");

  const auto builder = [segments](const Frame& frame, std::vector<double>& xx, std::vector<double>& yy, std::vector<double>& zz) {
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");
#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
    for (std::size_t k = 0; k < segments.size(); ++k) {
      const auto [i, j] = segments[k];
      const Vec3 u = safe_unit(Vec3{xu[j] - xu[i], yu[j] - yu[i], zu[j] - zu[i]});
      xx[k] = u.x; yy[k] = u.y; zz[k] = u.z;
    }
  };

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("segment_reorientation.dat"))).lexically_normal();

  OrientationCorrMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = resolve_exact_frame_end(parse_correlator_spec(cfg, section, env.dt), env.follow,
                                                opt.range.frame_start,
                                                cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1)),
                                                env.input_path);
  opt.range.dry_run = env.dry_run;
  opt.corr = parse_correlator_spec(cfg, section, env.dt);
  opt.legendre_order = static_cast<int>(cfg.get_int64(section, "legendre_order", std::optional<std::int64_t>(2)));
  opt.object_kind = "polymer_segment";
  return std::make_unique<OrientationCorrMeasure>("segment_reorientation", instance, out.string(), sel, segments.size(), opt, builder);
}

MeasureCapabilities ring_normal_acf_caps(const IniConfig& cfg,
                                         const std::string& section,
                                         const std::string& instance,
                                         const MeasureBuildEnv& env) {
  (void)instance; (void)env;
  MeasureCapabilities caps;
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_identity_consistent = true;
  caps.requires_dfields = {"xu", "yu", "zu",
                           cfg.get_string(section, "ring_id_field"),
                           cfg.get_string(section, "ring_pos_field")};
  caps.scale = ScaleCompatibility{true, true, true};
  caps.group_refs.push_back(cfg.get_string(section, "group", std::optional<std::string>("all")));
  return caps;
}

std::unique_ptr<IMeasure> ring_normal_acf_create(const IniConfig& cfg,
                                                 const std::string& section,
                                                 const std::string& instance,
                                                 const MeasureBuildEnv& env,
                                                 const SystemContext& sysctx) {
  (void)sysctx;
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("ring_normal_acf factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "ring_normal_acf");
  const auto ring_id = measure_ext::integer_like_field_to_i64(frame0, cfg.get_string(section, "ring_id_field"), true);
  const SelectedChains ring_groups = build_selected_chains(sel, ring_id);
  const auto ring_pos = frame0.require_dfield(cfg.get_string(section, "ring_pos_field"));
  const OrderedChains ordered = build_ordered_chains_from_chain_pos(ring_groups, ring_pos, true);
  std::vector<std::vector<std::size_t>> rings;
  for (const auto& ring : ordered.atom_by_chain) {
    if (ring.size() >= 3) rings.push_back(ring);
  }
  if (rings.empty()) throw std::runtime_error("ring_normal_acf: no rings with at least 3 ordered atoms");

  const auto builder = [rings](const Frame& frame, std::vector<double>& xx, std::vector<double>& yy, std::vector<double>& zz) {
    const auto xu = frame.require_dfield("xu");
    const auto yu = frame.require_dfield("yu");
    const auto zu = frame.require_dfield("zu");
#if PILOTS_HAS_OPENMP
#pragma omp parallel for
#endif
    for (std::size_t k = 0; k < rings.size(); ++k) {
      const auto& ring = rings[k];
      Vec3 normal{};
      for (std::size_t i = 0; i < ring.size(); ++i) {
        const std::size_t a = ring[i];
        const std::size_t b = ring[(i + 1) % ring.size()];
        const double xa = xu[a], ya = yu[a], za = zu[a];
        const double xb = xu[b], yb = yu[b], zb = zu[b];
        normal.x += (ya - yb) * (za + zb);
        normal.y += (za - zb) * (xa + xb);
        normal.z += (xa - xb) * (ya + yb);
      }
      normal = safe_unit(normal);
      xx[k] = normal.x; yy[k] = normal.y; zz[k] = normal.z;
    }
  };

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("ring_normal_acf.dat"))).lexically_normal();

  OrientationCorrMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = resolve_exact_frame_end(parse_correlator_spec(cfg, section, env.dt), env.follow,
                                                opt.range.frame_start,
                                                cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1)),
                                                env.input_path);
  opt.range.dry_run = env.dry_run;
  opt.corr = parse_correlator_spec(cfg, section, env.dt);
  opt.legendre_order = static_cast<int>(cfg.get_int64(section, "legendre_order", std::optional<std::int64_t>(2)));
  opt.object_kind = "ring_normal";
  return std::make_unique<OrientationCorrMeasure>("ring_normal_acf", instance, out.string(), sel, rings.size(), opt, builder);
}

MeasureCapabilities dipole_family_caps(const IniConfig& cfg,
                                       const std::string& section,
                                       const std::string& instance,
                                       const MeasureBuildEnv& env) {
  (void)instance;
  MeasureCapabilities caps;
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_identity_consistent = true;
  caps.requires_dfields = {"xu", "yu", "zu", cfg.get_string(section, "charge_field", std::optional<std::string>("q"))};
  caps.scale = ScaleCompatibility{true, true, true};
  caps.group_refs.push_back(cfg.get_string(section, "group", std::optional<std::string>("all")));
  if (cfg.has_key(section, "entity_id_field")) {
    caps.requires_dfields.push_back(cfg.get_string(section, "entity_id_field"));
  } else if (env.first_frame && env.first_frame->has_mol) {
    caps.requires_i64fields.push_back("mol");
  } else {
    caps.requires_topology_sections.push_back("bonds");
  }
  return caps;
}

std::unique_ptr<IMeasure> create_total_dipole_time(const IniConfig& cfg,
                                                   const std::string& section,
                                                   const std::string& instance,
                                                   const MeasureBuildEnv& env,
                                                   const SystemContext& sysctx,
                                                   std::string type_name,
                                                   bool write_spectrum) {
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error(type_name + " factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, type_name);
  const auto entity_id = entity_id_per_atom_from_config(cfg, section, frame0, sysctx);
  const SelectedChains entities = build_selected_chains(sel, entity_id);

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const std::string acf_name = cfg.get_string(section, "output_acf", std::optional<std::string>(type_name + std::string("_acf.dat")));
  const std::string spec_name = cfg.get_string(section, "output_spectrum", std::optional<std::string>(type_name + std::string("_spectrum.dat")));
  const fs::path out_acf = (output_dir / acf_name).lexically_normal();
  const fs::path out_spec = (output_dir / spec_name).lexically_normal();

  TotalDipoleTimeMeasure::Options opt;
  opt.range.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.range.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.range.dry_run = env.dry_run;
  opt.max_lag_frames = cfg.get_int64(section, "max_lag_frames", std::optional<std::int64_t>(-1));
  opt.lag_stride = cfg.get_int64(section, "lag_stride", std::optional<std::int64_t>(1));
  opt.origin_stride = cfg.get_int64(section, "origin_stride", std::optional<std::int64_t>(1));
  opt.charge_field = cfg.get_string(section, "charge_field", std::optional<std::string>("q"));
  opt.write_spectrum = write_spectrum;
  if (write_spectrum) {
    opt.omega_list = measure_ext::parse_double_list(cfg, section, "omega_list");
    opt.delta_epsilon = cfg.get_double(section, "delta_epsilon", std::optional<double>(std::numeric_limits<double>::quiet_NaN()));
    opt.epsilon_infty = cfg.get_double(section, "epsilon_infty", std::optional<double>(std::numeric_limits<double>::quiet_NaN()));
  }
  return std::make_unique<TotalDipoleTimeMeasure>(type_name, instance, out_acf.string(), out_spec.string(), sel, entities, opt);
}

std::unique_ptr<IMeasure> dipole_acf_create(const IniConfig& cfg,
                                            const std::string& section,
                                            const std::string& instance,
                                            const MeasureBuildEnv& env,
                                            const SystemContext& sysctx) {
  return create_total_dipole_time(cfg, section, instance, env, sysctx, "dipole_acf", false);
}

std::unique_ptr<IMeasure> dielectric_spectrum_create(const IniConfig& cfg,
                                                     const std::string& section,
                                                     const std::string& instance,
                                                     const MeasureBuildEnv& env,
                                                     const SystemContext& sysctx) {
  return create_total_dipole_time(cfg, section, instance, env, sysctx, "dielectric_spectrum", true);
}

static MeasureRegistrar g_reg_bond_vector_acf("bond_vector_acf", &bond_vector_acf_caps, &bond_vector_acf_create);
static MeasureRegistrar g_reg_segment_reorientation("segment_reorientation", &segment_reorientation_caps, &segment_reorientation_create);
static MeasureRegistrar g_reg_ring_normal_acf("ring_normal_acf", &ring_normal_acf_caps, &ring_normal_acf_create);
static MeasureRegistrar g_reg_dipole_acf("dipole_acf", &dipole_family_caps, &dipole_acf_create);
static MeasureRegistrar g_reg_dielectric_spectrum("dielectric_spectrum", &dipole_family_caps, &dielectric_spectrum_create);

} // namespace
} // namespace pilots
