#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "MeasureCommon.hpp"
#include "pilots/measures/IMeasure.hpp"
#include "pilots/measures/MeasureRegistry.hpp"
#include "pilots/select/SelectionView.hpp"
#include "pilots/util/AtomicFile.hpp"

namespace fs = std::filesystem;

namespace pilots {
namespace {

using measure_ext::Axis1D;
using measure_ext::SelectedChains;
using measure_ext::axis1d_name;
using measure_ext::axis_length;
using measure_ext::box_volume;
using measure_ext::build_selected_chains;
using measure_ext::dstr;
using measure_ext::entity_id_per_atom_from_config;
using measure_ext::get_static_combined_view;
using measure_ext::orth_area_for_axis;
using measure_ext::parse_axis1d;
using measure_ext::primary_axis_coord;
using measure_ext::resolve_path;

constexpr double kBoltzmann = 1.380649e-23;
constexpr double kEps0 = 8.8541878128e-12;
constexpr double kElementaryCharge = 1.602176634e-19;

struct DielectricCommonOptions {
  std::int64_t frame_start = 0;
  std::int64_t frame_end = -1;
  std::string charge_field = "q";
  std::string length_unit = "angstrom";
  std::string charge_unit = "e";
  double temperature = 0.0;
  double entity_charge_tolerance = 1e-6;
  bool allow_charged_entities = false;
  bool dry_run = false;
};

inline double length_to_m(std::string s) {
  for (auto& c : s) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
  if (s == "a" || s == "ang" || s == "angstrom" || s == "angstroms") return 1e-10;
  if (s == "nm" || s == "nanometer" || s == "nanometers") return 1e-9;
  if (s == "m" || s == "meter" || s == "meters") return 1.0;
  throw std::runtime_error("invalid length_unit: '" + s + "' (use angstrom|nm|m)");
}

inline double charge_to_c(std::string s) {
  for (auto& c : s) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
  if (s == "e" || s == "elem" || s == "elementary") return kElementaryCharge;
  if (s == "c" || s == "coulomb" || s == "coulombs") return 1.0;
  throw std::runtime_error("invalid charge_unit: '" + s + "' (use e|coulomb)");
}

inline bool frame_in_range(std::size_t frame_index,
                           std::int64_t frame_start,
                           std::int64_t frame_end) {
  const std::int64_t fi = static_cast<std::int64_t>(frame_index);
  if (fi < frame_start) return false;
  if (frame_end >= 0 && fi > frame_end) return false;
  return true;
}

inline double normal_component(Axis1D axis, double x, double y, double z) {
  switch (axis) {
    case Axis1D::X: return x;
    case Axis1D::Y: return y;
    case Axis1D::Z: return z;
  }
  return z;
}

inline double parallel_average_component(Axis1D axis, double x, double y, double z) {
  switch (axis) {
    case Axis1D::X: return 0.5 * (y + z);
    case Axis1D::Y: return 0.5 * (x + z);
    case Axis1D::Z: return 0.5 * (x + y);
  }
  return 0.5 * (x + y);
}

struct EntityDipoleFrame {
  double mx = 0.0;
  double my = 0.0;
  double mz = 0.0;
  double volume = 0.0;
  double axis_length = 0.0;
  std::vector<double> bin_mx;
  std::vector<double> bin_my;
  std::vector<double> bin_mz;
};

inline EntityDipoleFrame compute_entity_dipoles(const Frame& frame,
                                                const SelectionView& sel,
                                                const SelectedChains& entities,
                                                const DielectricCommonOptions& opt,
                                                bool with_bins,
                                                Axis1D axis,
                                                std::size_t n_bins) {
  const auto q = frame.require_dfield(opt.charge_field);
  const auto xu = frame.require_dfield("xu");
  const auto yu = frame.require_dfield("yu");
  const auto zu = frame.require_dfield("zu");

  EntityDipoleFrame out;
  out.volume = box_volume(frame.box);
  out.axis_length = measure_ext::axis_length(frame.box, axis);
  if (with_bins) {
    out.bin_mx.assign(n_bins, 0.0);
    out.bin_my.assign(n_bins, 0.0);
    out.bin_mz.assign(n_bins, 0.0);
  }

  for (std::size_t c = 0; c < entities.n_chains(); ++c) {
    const auto& atoms = entities.atom_by_chain[c];
    const double invn = 1.0 / static_cast<double>(atoms.size());
    double cx = 0.0, cy = 0.0, cz = 0.0;
    double qsum = 0.0;
    for (const std::size_t i : atoms) {
      cx += xu[i];
      cy += yu[i];
      cz += zu[i];
      qsum += q[i];
    }
    cx *= invn;
    cy *= invn;
    cz *= invn;
    if (!opt.allow_charged_entities && std::abs(qsum) > opt.entity_charge_tolerance) {
      throw std::runtime_error("dielectric measure requires approximately neutral grouped entities; entity id " +
                               std::to_string(entities.chain_id[c]) + " has net charge " + dstr(qsum));
    }

    double mux = 0.0, muy = 0.0, muz = 0.0;
    for (const std::size_t i : atoms) {
      mux += q[i] * (xu[i] - cx);
      muy += q[i] * (yu[i] - cy);
      muz += q[i] * (zu[i] - cz);
    }
    out.mx += mux;
    out.my += muy;
    out.mz += muz;

    if (with_bins) {
      const double s = primary_axis_coord(frame.box, cx, cy, cz, axis) / out.axis_length;
      std::size_t b = static_cast<std::size_t>(std::floor(s * static_cast<double>(n_bins)));
      if (b >= n_bins) b = n_bins - 1;
      out.bin_mx[b] += mux;
      out.bin_my[b] += muy;
      out.bin_mz[b] += muz;
    }
  }

  return out;
}

class BulkDielectricMeasure final : public IMeasure {
public:
  BulkDielectricMeasure(std::string instance_name,
                        std::string output_path,
                        SelectionView sel,
                        SelectedChains entities,
                        DielectricCommonOptions opt)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        entities_(std::move(entities)),
        opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    validate_common_();
    if (entities_.n_chains() == 0) throw std::runtime_error("bulk_dielectric: no grouped entities are available");
  }

  std::string type() const override { return "bulk_dielectric"; }
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
    od.columns = {"n_frames", "volume_mean", "mx_mean", "my_mean", "mz_mean", "var_mx", "var_my", "var_mz", "eps_x", "eps_y", "eps_z", "eps_iso"};
    md.outputs.push_back(std::move(od));

    md.params["temperature"] = dstr(opt_.temperature);
    md.params["charge_field"] = opt_.charge_field;
    md.params["length_unit"] = opt_.length_unit;
    md.params["charge_unit"] = opt_.charge_unit;
    md.params["frame_start"] = std::to_string(opt_.frame_start);
    md.params["frame_end"] = std::to_string(opt_.frame_end);
    md.params["n_entities"] = std::to_string(entities_.n_chains());
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    started_ = true;
    if (!opt_.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.frame_start, opt_.frame_end)) return;
    frames_.push_back(compute_entity_dipoles(frame, sel_, entities_, opt_, false, Axis1D::Z, 1));
  }

  void flush_partial() override {
    if (opt_.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      write_header_(ofs);
      if (frames_.empty()) {
        ofs << "0 0 0 0 0 0 0 0 0 0 0 0\n";
        return;
      }
      const double lscale = length_to_m(opt_.length_unit);
      const double qscale = charge_to_c(opt_.charge_unit);
      const double dscale = lscale * qscale;
      const double vscale = lscale * lscale * lscale;

      double sum_mx = 0.0, sum_my = 0.0, sum_mz = 0.0, sum_v = 0.0;
      double sum_mx2 = 0.0, sum_my2 = 0.0, sum_mz2 = 0.0;
      for (const auto& f : frames_) {
        const double mx = f.mx * dscale;
        const double my = f.my * dscale;
        const double mz = f.mz * dscale;
        const double v = f.volume * vscale;
        sum_mx += mx;
        sum_my += my;
        sum_mz += mz;
        sum_v += v;
        sum_mx2 += mx * mx;
        sum_my2 += my * my;
        sum_mz2 += mz * mz;
      }
      const double invn = 1.0 / static_cast<double>(frames_.size());
      const double mean_mx = sum_mx * invn;
      const double mean_my = sum_my * invn;
      const double mean_mz = sum_mz * invn;
      const double mean_v = sum_v * invn;
      const double var_mx = std::max(0.0, sum_mx2 * invn - mean_mx * mean_mx);
      const double var_my = std::max(0.0, sum_my2 * invn - mean_my * mean_my);
      const double var_mz = std::max(0.0, sum_mz2 * invn - mean_mz * mean_mz);
      const double pref = 1.0 / (kEps0 * mean_v * kBoltzmann * opt_.temperature);
      const double eps_x = 1.0 + pref * var_mx;
      const double eps_y = 1.0 + pref * var_my;
      const double eps_z = 1.0 + pref * var_mz;
      const double eps_iso = 1.0 + pref * (var_mx + var_my + var_mz) / 3.0;

      ofs << frames_.size() << ' '
          << std::setprecision(17) << mean_v << ' '
          << std::setprecision(17) << mean_mx << ' '
          << std::setprecision(17) << mean_my << ' '
          << std::setprecision(17) << mean_mz << ' '
          << std::setprecision(17) << var_mx << ' '
          << std::setprecision(17) << var_my << ' '
          << std::setprecision(17) << var_mz << ' '
          << std::setprecision(17) << eps_x << ' '
          << std::setprecision(17) << eps_y << ' '
          << std::setprecision(17) << eps_z << ' '
          << std::setprecision(17) << eps_iso << '\n';
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
  SelectionView sel_;
  SelectedChains entities_;
  DielectricCommonOptions opt_;
  std::vector<EntityDipoleFrame> frames_;
  bool started_ = false;

  void validate_common_() const {
    if (opt_.frame_start < 0) throw std::runtime_error("bulk_dielectric: frame_start must be >= 0");
    if (opt_.frame_end >= 0 && opt_.frame_end < opt_.frame_start) {
      throw std::runtime_error("bulk_dielectric: frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.temperature > 0.0)) throw std::runtime_error("bulk_dielectric: temperature must be > 0");
    if (sel_.idx.empty()) throw std::runtime_error("bulk_dielectric: selection is empty");
  }

  void write_header_(std::ostream& ofs) const {
    ofs << "# PILOTS: bulk dielectric constant from total neutral-entity dipole fluctuations\n";
    ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
    ofs << "# n_entities: " << entities_.n_chains() << "\n";
    ofs << "# charge_field: " << opt_.charge_field << "\n";
    ofs << "# charge_unit: " << opt_.charge_unit << ", length_unit: " << opt_.length_unit << "\n";
    ofs << "# temperature: " << std::setprecision(17) << opt_.temperature << " K\n";
    ofs << "# frame_range: [" << opt_.frame_start << ", " << opt_.frame_end << "]\n";
    ofs << "# columns: n_frames  volume_mean  mx_mean  my_mean  mz_mean  var_mx  var_my  var_mz  eps_x  eps_y  eps_z  eps_iso\n";
  }
};

class SlabDielectricMeasure final : public IMeasure {
public:
  struct Options {
    DielectricCommonOptions common;
    Axis1D axis = Axis1D::Z;
    std::size_t n_bins = 0;
  };

  SlabDielectricMeasure(std::string instance_name,
                        std::string output_path,
                        SelectionView sel,
                        SelectedChains entities,
                        Options opt)
      : instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        entities_(std::move(entities)),
        opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    sel_ = SelectionView{sel_name_owned_, sel.idx};
    validate_common_();
    if (entities_.n_chains() == 0) throw std::runtime_error("slab_dielectric: no grouped entities are available");
  }

  std::string type() const override { return "slab_dielectric"; }
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
    od.x_axis = "bin";
    od.x_unit = "fractional_slab_coordinate";
    od.columns = {"bin", "coord_lo_frac", "coord_hi_frac", "coord_center_frac", "coord_center_mean", "eps_parallel", "eps_perp", "chi_parallel", "chi_perp", "cov_parallel", "cov_perp", "n_frames"};
    md.outputs.push_back(std::move(od));

    md.params["axis"] = axis1d_name(opt_.axis);
    md.params["n_bins"] = std::to_string(opt_.n_bins);
    md.params["temperature"] = dstr(opt_.common.temperature);
    md.params["charge_field"] = opt_.common.charge_field;
    md.params["length_unit"] = opt_.common.length_unit;
    md.params["charge_unit"] = opt_.common.charge_unit;
    md.params["frame_start"] = std::to_string(opt_.common.frame_start);
    md.params["frame_end"] = std::to_string(opt_.common.frame_end);
    md.params["n_entities"] = std::to_string(entities_.n_chains());
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    started_ = true;
    if (!opt_.common.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.common.frame_start, opt_.common.frame_end)) return;
    frames_.push_back(compute_entity_dipoles(frame, sel_, entities_, opt_.common, true, opt_.axis, opt_.n_bins));
  }

  void flush_partial() override {
    if (opt_.common.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      write_header_(ofs);
      if (frames_.empty()) return;
      const double lscale = length_to_m(opt_.common.length_unit);
      const double qscale = charge_to_c(opt_.common.charge_unit);
      const double dscale = lscale * qscale;
      double mean_vbin = 0.0;
      double mean_axis_length = 0.0;
      double mean_mx = 0.0, mean_my = 0.0, mean_mz = 0.0;
      for (const auto& f : frames_) {
        mean_vbin += (f.volume / static_cast<double>(opt_.n_bins));
        mean_axis_length += f.axis_length;
        mean_mx += f.mx;
        mean_my += f.my;
        mean_mz += f.mz;
      }
      const double invn = 1.0 / static_cast<double>(frames_.size());
      mean_vbin *= invn * lscale * lscale * lscale;
      mean_axis_length *= invn * lscale;
      mean_mx *= invn * dscale;
      mean_my *= invn * dscale;
      mean_mz *= invn * dscale;

      for (std::size_t b = 0; b < opt_.n_bins; ++b) {
        double mean_bx = 0.0, mean_by = 0.0, mean_bz = 0.0;
        double mean_bxMx = 0.0, mean_byMy = 0.0, mean_bzMz = 0.0;
        for (const auto& f : frames_) {
          const double bx = f.bin_mx[b] * dscale;
          const double by = f.bin_my[b] * dscale;
          const double bz = f.bin_mz[b] * dscale;
          const double mx = f.mx * dscale;
          const double my = f.my * dscale;
          const double mz = f.mz * dscale;
          mean_bx += bx;
          mean_by += by;
          mean_bz += bz;
          mean_bxMx += bx * mx;
          mean_byMy += by * my;
          mean_bzMz += bz * mz;
        }
        mean_bx *= invn;
        mean_by *= invn;
        mean_bz *= invn;
        mean_bxMx *= invn;
        mean_byMy *= invn;
        mean_bzMz *= invn;

        const double cov_x = mean_bxMx - mean_bx * mean_mx;
        const double cov_y = mean_byMy - mean_by * mean_my;
        const double cov_z = mean_bzMz - mean_bz * mean_mz;
        const double cov_parallel = parallel_average_component(opt_.axis, cov_x, cov_y, cov_z);
        const double cov_perp = normal_component(opt_.axis, cov_x, cov_y, cov_z);
        const double chi_parallel = cov_parallel / (kEps0 * mean_vbin * kBoltzmann * opt_.common.temperature);
        const double chi_perp = cov_perp / (kEps0 * mean_vbin * kBoltzmann * opt_.common.temperature);
        const double eps_parallel = 1.0 + chi_parallel;
        const double denom = 1.0 - chi_perp;
        const double eps_perp = (std::abs(denom) > 1e-12) ? (1.0 / denom) : ((denom >= 0.0) ? 1e12 : -1e12);

        const double lo_frac = static_cast<double>(b) / static_cast<double>(opt_.n_bins);
        const double hi_frac = static_cast<double>(b + 1) / static_cast<double>(opt_.n_bins);
        const double center_frac = 0.5 * (lo_frac + hi_frac);
        const double center = center_frac * mean_axis_length;

        ofs << b << ' '
            << std::setprecision(17) << lo_frac << ' '
            << std::setprecision(17) << hi_frac << ' '
            << std::setprecision(17) << center_frac << ' '
            << std::setprecision(17) << center << ' '
            << std::setprecision(17) << eps_parallel << ' '
            << std::setprecision(17) << eps_perp << ' '
            << std::setprecision(17) << chi_parallel << ' '
            << std::setprecision(17) << chi_perp << ' '
            << std::setprecision(17) << cov_parallel << ' '
            << std::setprecision(17) << cov_perp << ' '
            << frames_.size() << '\n';
      }
    });
  }

  void finalize() override {
    if (opt_.common.dry_run) return;
    flush_partial();
  }

private:
  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  SelectionView sel_;
  SelectedChains entities_;
  Options opt_;
  std::vector<EntityDipoleFrame> frames_;
  bool started_ = false;

  void validate_common_() const {
    if (opt_.common.frame_start < 0) throw std::runtime_error("slab_dielectric: frame_start must be >= 0");
    if (opt_.common.frame_end >= 0 && opt_.common.frame_end < opt_.common.frame_start) {
      throw std::runtime_error("slab_dielectric: frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.common.temperature > 0.0)) throw std::runtime_error("slab_dielectric: temperature must be > 0");
    if (opt_.n_bins == 0) throw std::runtime_error("slab_dielectric: n_bins must be >= 1");
    if (sel_.idx.empty()) throw std::runtime_error("slab_dielectric: selection is empty");
  }

  void write_header_(std::ostream& ofs) const {
    ofs << "# PILOTS: slab/local dielectric profile from neutral-entity dipole fluctuations\n";
    ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
    ofs << "# n_entities: " << entities_.n_chains() << "\n";
    ofs << "# axis: " << axis1d_name(opt_.axis) << ", n_bins: " << opt_.n_bins << "\n";
    ofs << "# charge_field: " << opt_.common.charge_field << "\n";
    ofs << "# charge_unit: " << opt_.common.charge_unit << ", length_unit: " << opt_.common.length_unit << "\n";
    ofs << "# temperature: " << std::setprecision(17) << opt_.common.temperature << " K\n";
    ofs << "# frame_range: [" << opt_.common.frame_start << ", " << opt_.common.frame_end << "]\n";
    ofs << "# columns: bin  coord_lo_frac  coord_hi_frac  coord_center_frac  coord_center_mean  eps_parallel  eps_perp  chi_parallel  chi_perp  cov_parallel  cov_perp  n_frames\n";
  }
};

void append_dielectric_caps(const IniConfig& cfg,
                            const std::string& section,
                            const MeasureBuildEnv& env,
                            MeasureCapabilities& caps) {
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_identity_consistent = true;
  caps.requires_dfields = {"xu", "yu", "zu", cfg.get_string(section, "charge_field", std::optional<std::string>("q"))};
  caps.scale = ScaleCompatibility{true, true, true};
  caps.group_refs.push_back(cfg.get_string(section, "group", std::optional<std::string>("all")));
  if (cfg.has_key(section, "entity_id_field")) {
    caps.requires_dfields.push_back(cfg.get_string(section, "entity_id_field"));
    return;
  }
  if (env.first_frame && env.first_frame->has_mol) {
    caps.requires_i64fields.push_back("mol");
    return;
  }
  caps.requires_topology_sections.push_back("bonds");
}

MeasureCapabilities bulk_dielectric_caps(const IniConfig& cfg,
                                         const std::string& section,
                                         const std::string& instance,
                                         const MeasureBuildEnv& env) {
  (void)instance;
  MeasureCapabilities caps;
  append_dielectric_caps(cfg, section, env, caps);
  return caps;
}

std::unique_ptr<IMeasure> bulk_dielectric_create(const IniConfig& cfg,
                                                 const std::string& section,
                                                 const std::string& instance,
                                                 const MeasureBuildEnv& env,
                                                 const SystemContext& sysctx) {
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("bulk_dielectric factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "bulk_dielectric");
  const auto entity_ids = entity_id_per_atom_from_config(cfg, section, frame0, sysctx);
  SelectedChains entities = build_selected_chains(sel, entity_ids);

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("bulk_dielectric.dat"))).lexically_normal();

  DielectricCommonOptions opt;
  opt.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.charge_field = cfg.get_string(section, "charge_field", std::optional<std::string>("q"));
  opt.length_unit = cfg.get_string(section, "length_unit", std::optional<std::string>("angstrom"));
  opt.charge_unit = cfg.get_string(section, "charge_unit", std::optional<std::string>("e"));
  opt.temperature = cfg.get_double(section, "temperature");
  opt.entity_charge_tolerance = cfg.get_double(section, "entity_charge_tolerance", std::optional<double>(1e-6));
  opt.allow_charged_entities = cfg.get_bool(section, "allow_charged_entities", std::optional<bool>(false));
  opt.dry_run = env.dry_run;

  return std::make_unique<BulkDielectricMeasure>(instance, out.string(), sel, std::move(entities), opt);
}

MeasureCapabilities slab_dielectric_caps(const IniConfig& cfg,
                                         const std::string& section,
                                         const std::string& instance,
                                         const MeasureBuildEnv& env) {
  (void)instance;
  MeasureCapabilities caps;
  append_dielectric_caps(cfg, section, env, caps);
  return caps;
}

std::unique_ptr<IMeasure> slab_dielectric_create(const IniConfig& cfg,
                                                 const std::string& section,
                                                 const std::string& instance,
                                                 const MeasureBuildEnv& env,
                                                 const SystemContext& sysctx) {
  if (!env.first_frame || !env.selection_provider) throw std::runtime_error("slab_dielectric factory: missing first_frame or SelectionProvider");
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, "slab_dielectric");
  const auto entity_ids = entity_id_per_atom_from_config(cfg, section, frame0, sysctx);
  SelectedChains entities = build_selected_chains(sel, entity_ids);

  const fs::path output_dir = cfg.get_string(section, "output_dir", std::optional<std::string>("")).empty()
      ? env.output_dir_general
      : resolve_path(env.cfg_dir, cfg.get_string(section, "output_dir"));
  if (!env.dry_run) fs::create_directories(output_dir);
  const fs::path out = (output_dir / cfg.get_string(section, "output", std::optional<std::string>("slab_dielectric.dat"))).lexically_normal();

  SlabDielectricMeasure::Options opt;
  opt.common.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.common.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.common.charge_field = cfg.get_string(section, "charge_field", std::optional<std::string>("q"));
  opt.common.length_unit = cfg.get_string(section, "length_unit", std::optional<std::string>("angstrom"));
  opt.common.charge_unit = cfg.get_string(section, "charge_unit", std::optional<std::string>("e"));
  opt.common.temperature = cfg.get_double(section, "temperature");
  opt.common.entity_charge_tolerance = cfg.get_double(section, "entity_charge_tolerance", std::optional<double>(1e-6));
  opt.common.allow_charged_entities = cfg.get_bool(section, "allow_charged_entities", std::optional<bool>(false));
  opt.common.dry_run = env.dry_run;
  opt.axis = parse_axis1d(cfg.get_string(section, "axis", std::optional<std::string>("z")));
  opt.n_bins = static_cast<std::size_t>(cfg.get_int64(section, "n_bins", std::optional<std::int64_t>(100)));

  return std::make_unique<SlabDielectricMeasure>(instance, out.string(), sel, std::move(entities), opt);
}

static MeasureRegistrar g_register_bulk_dielectric("bulk_dielectric", &bulk_dielectric_caps, &bulk_dielectric_create);
static MeasureRegistrar g_register_slab_dielectric("slab_dielectric", &slab_dielectric_caps, &slab_dielectric_create);
static MeasureRegistrar g_register_slab_dielectric_alias("layered_dielectric", &slab_dielectric_caps, &slab_dielectric_create);

} // namespace
} // namespace pilots
