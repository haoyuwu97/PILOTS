#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <limits>
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
using measure_ext::SlabFrameTransform;
using measure_ext::SlabTransformOptions;
using measure_ext::append_integer_like_field_cap;
using measure_ext::axis1d_name;
using measure_ext::build_optional_anchor_selection;
using measure_ext::build_selected_chains;
using measure_ext::box_volume;
using measure_ext::dstr;
using measure_ext::entity_id_per_atom_from_config;
using measure_ext::get_static_combined_view;
using measure_ext::get_static_group_view;
using measure_ext::make_slab_frame_transform;
using measure_ext::mass_by_atom_from_config;
using measure_ext::parse_axis1d;
using measure_ext::resolve_measure_output_path;
using measure_ext::slab_transform_requested;

constexpr double kBoltzmann = 1.380649e-23;
constexpr double kEps0 = 8.8541878128e-12;
constexpr double kElementaryCharge = 1.602176634e-19;
constexpr double kContributionFractionTolerance = 1e-12;

enum class ChargeModel {
  NeutralOnly,
  FragmentInternal,
};

enum class EntityReference {
  Geometric,
  Mass,
};

struct DielectricCommonOptions {
  std::int64_t frame_start = 0;
  std::int64_t frame_end = -1;
  std::string charge_field = "q";
  std::string length_unit = "angstrom";
  std::string charge_unit = "e";
  double temperature = 0.0;
  double entity_charge_tolerance = 1e-6;
  ChargeModel charge_model = ChargeModel::NeutralOnly;
  EntityReference entity_reference = EntityReference::Geometric;
  std::vector<double> mass_by_atom;
  bool dry_run = false;
};

inline std::string lower_copy(std::string s) {
  for (auto& c : s) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
  return s;
}

inline double length_to_m(std::string s) {
  s = lower_copy(std::move(s));
  if (s == "a" || s == "ang" || s == "angstrom" || s == "angstroms") return 1e-10;
  if (s == "nm" || s == "nanometer" || s == "nanometers") return 1e-9;
  if (s == "m" || s == "meter" || s == "meters") return 1.0;
  throw std::runtime_error("invalid length_unit: '" + s + "' (use angstrom|nm|m)");
}

inline double charge_to_c(std::string s) {
  s = lower_copy(std::move(s));
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

inline EntityReference parse_entity_reference(const std::string& raw) {
  const std::string s = lower_copy(raw);
  if (s == "geometric" || s == "geom" || s == "center") return EntityReference::Geometric;
  if (s == "mass" || s == "com" || s == "center_of_mass") return EntityReference::Mass;
  throw std::runtime_error("invalid entity_reference: '" + raw + "' (use geometric|mass)");
}

inline std::string entity_reference_name(EntityReference ref) {
  switch (ref) {
    case EntityReference::Geometric: return "geometric";
    case EntityReference::Mass: return "mass";
  }
  return "geometric";
}

inline std::string charge_model_name(ChargeModel mode) {
  switch (mode) {
    case ChargeModel::NeutralOnly: return "neutral_only";
    case ChargeModel::FragmentInternal: return "fragment_internal";
  }
  return "neutral_only";
}

inline ChargeModel parse_charge_model_impl(const std::string& raw) {
  const std::string s = lower_copy(raw);
  if (s == "neutral_only" || s == "neutral" || s == "neutral_entities") {
    return ChargeModel::NeutralOnly;
  }
  if (s == "fragment_internal" || s == "internal_fragment" || s == "charged_fragment"
      || s == "mixed_fragment" || s == "fragment") {
    return ChargeModel::FragmentInternal;
  }
  throw std::runtime_error(
      "invalid charge_model: '" + raw + "' (use neutral_only|fragment_internal)");
}

inline ChargeModel charge_model_from_config(const IniConfig& cfg,
                                            const std::string& section,
                                            ChargeModel default_mode) {
  if (cfg.has_key(section, "charge_model")) {
    return parse_charge_model_impl(cfg.get_string(section, "charge_model"));
  }
  if (cfg.has_key(section, "allow_charged_entities")) {
    return cfg.get_bool(section, "allow_charged_entities", std::optional<bool>(false))
               ? ChargeModel::FragmentInternal
               : ChargeModel::NeutralOnly;
  }
  return default_mode;
}

struct EntityDipoleFrame {
  double mx = 0.0;
  double my = 0.0;
  double mz = 0.0;
  double volume = 0.0;
  double domain_length = 0.0;
  std::vector<double> bin_mx;
  std::vector<double> bin_my;
  std::vector<double> bin_mz;
  double charged_entities = 0.0;
  double sum_abs_entity_charge = 0.0;
};

inline EntityDipoleFrame compute_entity_internal_dipoles(const Frame& frame,
                                                         const SelectedChains& entities,
                                                         const DielectricCommonOptions& opt,
                                                         const SlabFrameTransform* slab_tf,
                                                         Axis1D axis,
                                                         std::size_t n_bins) {
  const auto q = frame.require_dfield(opt.charge_field);
  const auto xu = frame.require_dfield("xu");
  const auto yu = frame.require_dfield("yu");
  const auto zu = frame.require_dfield("zu");

  EntityDipoleFrame out;
  out.volume = box_volume(frame.box);
  out.domain_length = slab_tf ? slab_tf->domain_length() : measure_ext::axis_length(frame.box, axis);
  if (n_bins > 0) {
    out.bin_mx.assign(n_bins, 0.0);
    out.bin_my.assign(n_bins, 0.0);
    out.bin_mz.assign(n_bins, 0.0);
  }

  for (std::size_t c = 0; c < entities.n_chains(); ++c) {
    const auto& atoms = entities.atom_by_chain[c];
    if (atoms.empty()) continue;

    double wsum = 0.0;
    double rx = 0.0, ry = 0.0, rz = 0.0;
    double qsum = 0.0;
    for (const std::size_t i : atoms) {
      const double w = (opt.entity_reference == EntityReference::Mass)
                         ? opt.mass_by_atom.at(i)
                         : 1.0;
      if (!(w > 0.0)) {
        throw std::runtime_error("dielectric measure requires positive entity weights; check mass_field/topology masses");
      }
      wsum += w;
      rx += w * xu[i];
      ry += w * yu[i];
      rz += w * zu[i];
      qsum += q[i];
    }
    if (!(wsum > 0.0)) {
      throw std::runtime_error("dielectric measure encountered an entity with non-positive total weight");
    }

    const double cx = rx / wsum;
    const double cy = ry / wsum;
    const double cz = rz / wsum;
    const double abs_qsum = std::abs(qsum);
    if (opt.charge_model == ChargeModel::NeutralOnly && abs_qsum > opt.entity_charge_tolerance) {
      throw std::runtime_error(
          "dielectric measure requires approximately neutral grouped entities in neutral_only mode; entity id "
          + std::to_string(entities.chain_id[c]) + " has net charge " + dstr(qsum)
          + ". Use charge_model=fragment_internal or slab_fragment_dielectric for charged mixtures.");
    }

    if (abs_qsum > opt.entity_charge_tolerance) {
      out.charged_entities += 1.0;
      out.sum_abs_entity_charge += abs_qsum;
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

    if (n_bins > 0) {
      double s = 0.0;
      if (slab_tf) {
        s = slab_tf->point_to_domain_fraction(frame.box, cx, cy, cz);
      } else {
        const double L = measure_ext::axis_length(frame.box, axis);
        s = measure_ext::primary_axis_coord(frame.box, cx, cy, cz, axis) / L;
      }
      s = std::clamp(s, 0.0, 1.0);
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
  BulkDielectricMeasure(std::string type_name,
                        std::string instance_name,
                        std::string output_path,
                        SelectionView sel,
                        SelectedChains entities,
                        DielectricCommonOptions opt)
      : type_name_(std::move(type_name)),
        instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        entities_(std::move(entities)),
        opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    sel_idx_owned_.assign(sel.idx.begin(), sel.idx.end());
    sel_ = SelectionView{sel_name_owned_, std::span<const std::size_t>(sel_idx_owned_.data(), sel_idx_owned_.size())};
    validate_common_();
    if (entities_.n_chains() == 0) throw std::runtime_error(type_name_ + ": no grouped entities are available");
    if (opt_.entity_reference == EntityReference::Mass && opt_.mass_by_atom.size() != sel_idx_owned_.size() && !opt_.mass_by_atom.empty()) {
      // mass_by_atom is stored per atom in frame, not per selection; size is validated in create.
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
    md.params["charge_model"] = charge_model_name(opt_.charge_model);
    md.params["entity_reference"] = entity_reference_name(opt_.entity_reference);
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    started_ = true;
    if (!opt_.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.frame_start, opt_.frame_end)) return;
    frames_.push_back(compute_entity_internal_dipoles(frame, entities_, opt_, nullptr, Axis1D::Z, 0));
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
      double sum_charged = 0.0, sum_abs_q = 0.0;
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
        sum_charged += f.charged_entities;
        sum_abs_q += f.sum_abs_entity_charge;
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

      ofs << "# mean_charged_entities: " << std::setprecision(17) << (sum_charged * invn) << "\n";
      ofs << "# mean_abs_entity_charge(" << opt_.charge_unit << "): " << std::setprecision(17) << (sum_abs_q * invn) << "\n";
      if (opt_.charge_model == ChargeModel::FragmentInternal) {
        ofs << "# note: fragment_internal mode computes a finite, conduction-blocked response from entity-internal dipoles; translational free-charge response is excluded by construction.\n";
      }

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
  std::string type_name_;
  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  std::vector<std::size_t> sel_idx_owned_;
  SelectionView sel_{{}, {}};
  SelectedChains entities_;
  DielectricCommonOptions opt_;
  std::vector<EntityDipoleFrame> frames_;
  bool started_ = false;

  void validate_common_() const {
    if (opt_.frame_start < 0) throw std::runtime_error(type_name_ + ": frame_start must be >= 0");
    if (opt_.frame_end >= 0 && opt_.frame_end < opt_.frame_start) {
      throw std::runtime_error(type_name_ + ": frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.temperature > 0.0)) throw std::runtime_error(type_name_ + ": temperature must be > 0");
    if (sel_.idx.empty()) throw std::runtime_error(type_name_ + ": selection is empty");
    if (opt_.entity_reference == EntityReference::Mass && opt_.mass_by_atom.empty()) {
      throw std::runtime_error(type_name_ + ": entity_reference=mass requires mass_field or topology masses");
    }
  }

  void write_header_(std::ostream& ofs) const {
    ofs << "# PILOTS: bulk dielectric from entity-internal dipole fluctuations\n";
    ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
    ofs << "# n_entities: " << entities_.n_chains() << "\n";
    ofs << "# charge_field: " << opt_.charge_field << "\n";
    ofs << "# charge_unit: " << opt_.charge_unit << ", length_unit: " << opt_.length_unit << "\n";
    ofs << "# temperature: " << std::setprecision(17) << opt_.temperature << " K\n";
    ofs << "# frame_range: [" << opt_.frame_start << ", " << opt_.frame_end << "]\n";
    ofs << "# charge_model: " << charge_model_name(opt_.charge_model) << "\n";
    ofs << "# entity_reference: " << entity_reference_name(opt_.entity_reference) << "\n";
    if (opt_.charge_model == ChargeModel::NeutralOnly) {
      ofs << "# meaning: classic neutral-entity dipole-fluctuation dielectric.\n";
    } else {
      ofs << "# meaning: charged-mixture fragment-internal dielectric; finite and origin-invariant at the entity level, excluding translational ionic conductivity.\n";
    }
    ofs << "# columns: n_frames  volume_mean  mx_mean  my_mean  mz_mean  var_mx  var_my  var_mz  eps_x  eps_y  eps_z  eps_iso\n";
  }
};

class SlabDielectricMeasure final : public IMeasure {
public:
  struct Options {
    DielectricCommonOptions common;
    Axis1D axis = Axis1D::Z;
    std::size_t n_bins = 0;
    SlabTransformOptions slab;
  };

  SlabDielectricMeasure(std::string type_name,
                        std::string instance_name,
                        std::string output_path,
                        SelectionView sel,
                        std::optional<SelectionView> anchor_sel,
                        SelectedChains entities,
                        Options opt)
      : type_name_(std::move(type_name)),
        instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        entities_(std::move(entities)),
        opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    sel_idx_owned_.assign(sel.idx.begin(), sel.idx.end());
    sel_ = SelectionView{sel_name_owned_, std::span<const std::size_t>(sel_idx_owned_.data(), sel_idx_owned_.size())};
    if (anchor_sel.has_value()) {
      anchor_name_owned_ = std::string(anchor_sel->name);
      anchor_idx_owned_.assign(anchor_sel->idx.begin(), anchor_sel->idx.end());
      anchor_sel_ = SelectionView{anchor_name_owned_, std::span<const std::size_t>(anchor_idx_owned_.data(), anchor_idx_owned_.size())};
    }
    validate_common_();
    if (entities_.n_chains() == 0) throw std::runtime_error(type_name_ + ": no grouped entities are available");
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
    od.x_axis = "bin";
    od.x_unit = "fractional_slab_coordinate";
    od.columns = {"bin", "coord_lo_frac", "coord_hi_frac", "coord_center_frac", "coord_center_mean", "eps_parallel", "eps_perp", "chi_parallel", "chi_perp", "cov_parallel", "cov_perp", "n_frames"};
    md.outputs.push_back(std::move(od));

    md.params["temperature"] = dstr(opt_.common.temperature);
    md.params["charge_field"] = opt_.common.charge_field;
    md.params["length_unit"] = opt_.common.length_unit;
    md.params["charge_unit"] = opt_.common.charge_unit;
    md.params["frame_start"] = std::to_string(opt_.common.frame_start);
    md.params["frame_end"] = std::to_string(opt_.common.frame_end);
    md.params["axis"] = axis1d_name(opt_.axis);
    md.params["n_bins"] = std::to_string(opt_.n_bins);
    md.params["n_entities"] = std::to_string(entities_.n_chains());
    md.params["charge_model"] = charge_model_name(opt_.common.charge_model);
    md.params["entity_reference"] = entity_reference_name(opt_.common.entity_reference);
    md.params["slab_align_recenter"] = opt_.slab.slab_align_recenter ? "true" : "false";
    md.params["halfcell_fold"] = opt_.slab.halfcell_fold ? "true" : "false";
    md.params["target_center_frac"] = dstr(opt_.slab.target_center_frac);
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    started_ = true;
    if (!opt_.common.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.common.frame_start, opt_.common.frame_end)) return;
    const SlabFrameTransform tf = make_slab_frame_transform(frame, anchor_sel_ ? &(*anchor_sel_) : nullptr, opt_.slab);
    frames_.push_back(compute_entity_internal_dipoles(frame, entities_, opt_.common, &tf, opt_.axis, opt_.n_bins));
  }

  void flush_partial() override {
    if (opt_.common.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      write_header_(ofs);
      if (frames_.empty()) return;
      const double lscale = length_to_m(opt_.common.length_unit);
      const double qscale = charge_to_c(opt_.common.charge_unit);
      const double dscale = lscale * qscale;
      const double vscale = lscale * lscale * lscale;
      const double frames_inv = 1.0 / static_cast<double>(frames_.size());

      double mean_domain_length = 0.0;
      double mean_charged = 0.0;
      double mean_abs_q = 0.0;
      for (const auto& f : frames_) {
        mean_domain_length += f.domain_length;
        mean_charged += f.charged_entities;
        mean_abs_q += f.sum_abs_entity_charge;
      }
      mean_domain_length *= frames_inv * lscale;
      mean_charged *= frames_inv;
      mean_abs_q *= frames_inv;

      ofs << "# mean_charged_entities: " << std::setprecision(17) << mean_charged << "\n";
      ofs << "# mean_abs_entity_charge(" << opt_.common.charge_unit << "): " << std::setprecision(17) << mean_abs_q << "\n";
      if (opt_.common.charge_model == ChargeModel::FragmentInternal) {
        ofs << "# note: fragment_internal mode is a finite, conduction-blocked layered dielectric constructed from entity-internal dipoles. It is suitable for charged mixtures when a molecular/fragment entity decomposition is meaningful.\n";
      }

      for (std::size_t b = 0; b < opt_.n_bins; ++b) {
        double sum_bin_parallel = 0.0;
        double sum_bin_perp = 0.0;
        double sum_bin_parallel_total = 0.0;
        double sum_bin_perp_total = 0.0;
        double sum_total_parallel = 0.0;
        double sum_total_perp = 0.0;
        double sum_vbin = 0.0;

        for (const auto& f : frames_) {
          const double mx = f.mx * dscale;
          const double my = f.my * dscale;
          const double mz = f.mz * dscale;
          const double bin_mx = f.bin_mx[b] * dscale;
          const double bin_my = f.bin_my[b] * dscale;
          const double bin_mz = f.bin_mz[b] * dscale;
          const double total_parallel = parallel_average_component(opt_.axis, mx, my, mz);
          const double total_perp = normal_component(opt_.axis, mx, my, mz);
          const double bin_parallel = parallel_average_component(opt_.axis, bin_mx, bin_my, bin_mz);
          const double bin_perp = normal_component(opt_.axis, bin_mx, bin_my, bin_mz);
          const double vbin = (f.volume * vscale) / static_cast<double>(opt_.n_bins);

          sum_total_parallel += total_parallel;
          sum_total_perp += total_perp;
          sum_bin_parallel += bin_parallel;
          sum_bin_perp += bin_perp;
          sum_bin_parallel_total += bin_parallel * total_parallel;
          sum_bin_perp_total += bin_perp * total_perp;
          sum_vbin += vbin;
        }

        const double mean_total_parallel = sum_total_parallel * frames_inv;
        const double mean_total_perp = sum_total_perp * frames_inv;
        const double mean_bin_parallel = sum_bin_parallel * frames_inv;
        const double mean_bin_perp = sum_bin_perp * frames_inv;
        const double cov_parallel = sum_bin_parallel_total * frames_inv - mean_bin_parallel * mean_total_parallel;
        const double cov_perp = sum_bin_perp_total * frames_inv - mean_bin_perp * mean_total_perp;
        const double mean_vbin = sum_vbin * frames_inv;

        const double chi_parallel = cov_parallel / (kEps0 * mean_vbin * kBoltzmann * opt_.common.temperature);
        const double chi_perp = cov_perp / (kEps0 * mean_vbin * kBoltzmann * opt_.common.temperature);
        const double eps_parallel = 1.0 + chi_parallel;
        const double denom = 1.0 - chi_perp;
        const double eps_perp = (std::abs(denom) > 1e-12)
                                  ? (1.0 / denom)
                                  : ((denom >= 0.0) ? 1e12 : -1e12);

        const double lo_frac = static_cast<double>(b) / static_cast<double>(opt_.n_bins);
        const double hi_frac = static_cast<double>(b + 1) / static_cast<double>(opt_.n_bins);
        const double center_frac = 0.5 * (lo_frac + hi_frac);
        const double center = center_frac * mean_domain_length;

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
  std::string type_name_;
  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  std::vector<std::size_t> sel_idx_owned_;
  SelectionView sel_{{}, {}};
  std::string anchor_name_owned_;
  std::vector<std::size_t> anchor_idx_owned_;
  std::optional<SelectionView> anchor_sel_;
  SelectedChains entities_;
  Options opt_;
  std::vector<EntityDipoleFrame> frames_;
  bool started_ = false;

  void validate_common_() const {
    if (opt_.common.frame_start < 0) throw std::runtime_error(type_name_ + ": frame_start must be >= 0");
    if (opt_.common.frame_end >= 0 && opt_.common.frame_end < opt_.common.frame_start) {
      throw std::runtime_error(type_name_ + ": frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.common.temperature > 0.0)) throw std::runtime_error(type_name_ + ": temperature must be > 0");
    if (opt_.n_bins == 0) throw std::runtime_error(type_name_ + ": n_bins must be >= 1");
    if (sel_.idx.empty()) throw std::runtime_error(type_name_ + ": selection is empty");
    if (opt_.common.entity_reference == EntityReference::Mass && opt_.common.mass_by_atom.empty()) {
      throw std::runtime_error(type_name_ + ": entity_reference=mass requires mass_field or topology masses");
    }
  }

  void write_header_(std::ostream& ofs) const {
    ofs << "# PILOTS: slab dielectric profile from entity-internal dipole fluctuations\n";
    ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
    ofs << "# n_entities: " << entities_.n_chains() << "\n";
    ofs << "# axis: " << axis1d_name(opt_.axis) << ", n_bins: " << opt_.n_bins << "\n";
    ofs << "# charge_field: " << opt_.common.charge_field << "\n";
    ofs << "# charge_unit: " << opt_.common.charge_unit << ", length_unit: " << opt_.common.length_unit << "\n";
    ofs << "# temperature: " << std::setprecision(17) << opt_.common.temperature << " K\n";
    ofs << "# frame_range: [" << opt_.common.frame_start << ", " << opt_.common.frame_end << "]\n";
    ofs << "# charge_model: " << charge_model_name(opt_.common.charge_model) << "\n";
    ofs << "# entity_reference: " << entity_reference_name(opt_.common.entity_reference) << "\n";
    ofs << "# slab_align_recenter: " << (opt_.slab.slab_align_recenter ? "true" : "false")
        << ", halfcell_fold: " << (opt_.slab.halfcell_fold ? "true" : "false")
        << ", target_center_frac: " << std::setprecision(17) << opt_.slab.target_center_frac << "\n";
    if (anchor_sel_) {
      ofs << "# anchor_selection: " << anchor_sel_->name << " (n=" << anchor_sel_->idx.size() << ")\n";
    }
    if (opt_.common.charge_model == ChargeModel::NeutralOnly) {
      ofs << "# meaning: classic neutral-entity local dielectric profile.\n";
    } else {
      ofs << "# meaning: charged-mixture fragment-internal local dielectric profile; finite and conduction-blocked, excluding translational free-charge response.\n";
    }
    ofs << "# columns: bin  coord_lo_frac  coord_hi_frac  coord_center_frac  coord_center_mean  eps_parallel  eps_perp  chi_parallel  chi_perp  cov_parallel  cov_perp  n_frames\n";
  }
};




void append_dielectric_caps(const IniConfig& cfg,
                            const std::string& section,
                            const MeasureBuildEnv& env,
                            MeasureCapabilities& caps);
DielectricCommonOptions common_options_from_config(const IniConfig& cfg,
                                                    const std::string& section,
                                                    const Frame& frame0,
                                                    const SystemContext& sysctx,
                                                    ChargeModel default_model,
                                                    bool dry_run);
struct SpeciesDecomposition {
  std::vector<std::string> labels;
  std::vector<std::size_t> selected_counts;
  std::vector<int> atom_species;
  std::vector<std::uint8_t> in_background;
  std::size_t assigned_atoms = 0;
  std::size_t unassigned_atoms = 0;
  bool allow_partial_cover = false;
};

struct BackgroundDielectricFrame {
  double volume = 0.0;
  double domain_length = 0.0;
  std::vector<double> total_mx, total_my, total_mz;
  std::vector<std::vector<double>> bin_mx, bin_my, bin_mz;
};

inline std::vector<std::string> get_required_name_list(const IniConfig& cfg,
                                                        const std::string& section,
                                                        const std::string& key_name) {
  const auto out = cfg.get_list(section, key_name);
  if (out.empty()) {
    throw std::runtime_error("dielectric measure requires a non-empty '" + key_name + "' list");
  }
  return out;
}

inline double contribution_fraction_or_nan(double chi_contrib,
                                           double chi_total) {
  return (std::abs(chi_total) > kContributionFractionTolerance)
           ? (chi_contrib / chi_total)
           : std::numeric_limits<double>::quiet_NaN();
}

inline SpeciesDecomposition build_species_decomposition(const IniConfig& cfg,
                                                        const std::string& section,
                                                        const Frame& frame0,
                                                        SelectionProvider& sp,
                                                        const SelectionView& base_sel,
                                                        const std::string& who) {
  if (!cfg.has_key(section, "species_groups")) {
    throw std::runtime_error(who + ": species_groups is required");
  }
  const auto labels = get_required_name_list(cfg, section, "species_groups");
  const auto bg_labels = cfg.has_key(section, "background_groups")
                           ? get_required_name_list(cfg, section, "background_groups")
                           : labels;

  SpeciesDecomposition out;
  out.labels = labels;
  out.selected_counts.assign(labels.size(), 0);
  out.atom_species.assign(frame0.natoms, -1);
  out.in_background.assign(labels.size(), 0);
  out.allow_partial_cover = cfg.get_bool(section, "allow_partial_species_cover", std::optional<bool>(false));

  std::vector<std::uint8_t> in_base(frame0.natoms, 0);
  for (const std::size_t atom : base_sel.idx) {
    if (atom >= in_base.size()) throw std::runtime_error(who + ": selected atom index out of range");
    in_base[atom] = 1;
  }

  for (std::size_t s = 0; s < labels.size(); ++s) {
    const SelectionView gsel = get_static_group_view(sp, frame0, labels[s], who);
    std::size_t count = 0;
    for (const std::size_t atom : gsel.idx) {
      if (atom >= in_base.size()) throw std::runtime_error(who + ": group atom index out of range");
      if (!in_base[atom]) continue;
      if (out.atom_species[atom] >= 0 && out.atom_species[atom] != static_cast<int>(s)) {
        throw std::runtime_error(who + ": species_groups must be disjoint inside the base selection; atom "
                                 + std::to_string(atom) + " appears in multiple species groups");
      }
      if (out.atom_species[atom] != static_cast<int>(s)) {
        out.atom_species[atom] = static_cast<int>(s);
        ++count;
      }
    }
    if (count == 0) {
      throw std::runtime_error(who + ": species group '" + labels[s]
                               + "' has no atoms inside the base selection");
    }
    out.selected_counts[s] = count;
    out.assigned_atoms += count;
  }

  for (const std::size_t atom : base_sel.idx) {
    if (out.atom_species[atom] < 0) ++out.unassigned_atoms;
  }
  if (out.unassigned_atoms > 0 && !out.allow_partial_cover) {
    throw std::runtime_error(
        who + ": species_groups do not cover the full base selection; "
        + std::to_string(out.unassigned_atoms)
        + " selected atoms are unassigned. Add the missing atoms to species_groups or set allow_partial_species_cover=true to keep a partial decomposition.");
  }

  for (const auto& name : bg_labels) {
    bool found = false;
    for (std::size_t s = 0; s < labels.size(); ++s) {
      if (labels[s] == name) {
        out.in_background[s] = 1;
        found = true;
        break;
      }
    }
    if (!found) {
      throw std::runtime_error(who + ": background_groups entry '" + name
                               + "' is not present in species_groups");
    }
  }

  bool any_bg = false;
  for (const auto flag : out.in_background) {
    if (flag) { any_bg = true; break; }
  }
  if (!any_bg) throw std::runtime_error(who + ": background_groups selects no species");

  return out;
}

inline BackgroundDielectricFrame compute_background_dielectric_frame(
    const Frame& frame,
    const SelectedChains& entities,
    const DielectricCommonOptions& opt,
    const SlabFrameTransform* slab_tf,
    Axis1D axis,
    std::size_t n_bins,
    const SpeciesDecomposition& species) {
  const auto q = frame.require_dfield(opt.charge_field);
  const auto xu = frame.require_dfield("xu");
  const auto yu = frame.require_dfield("yu");
  const auto zu = frame.require_dfield("zu");

  BackgroundDielectricFrame out;
  out.volume = box_volume(frame.box);
  out.domain_length = slab_tf ? slab_tf->domain_length() : measure_ext::axis_length(frame.box, axis);
  const std::size_t ns = species.labels.size();
  out.total_mx.assign(ns, 0.0);
  out.total_my.assign(ns, 0.0);
  out.total_mz.assign(ns, 0.0);
  out.bin_mx.assign(ns, std::vector<double>(n_bins, 0.0));
  out.bin_my.assign(ns, std::vector<double>(n_bins, 0.0));
  out.bin_mz.assign(ns, std::vector<double>(n_bins, 0.0));

  for (std::size_t c = 0; c < entities.n_chains(); ++c) {
    const auto& atoms = entities.atom_by_chain[c];
    if (atoms.empty()) continue;

    double wsum = 0.0;
    double rx = 0.0, ry = 0.0, rz = 0.0;
    double qsum = 0.0;
    for (const std::size_t i : atoms) {
      const double w = (opt.entity_reference == EntityReference::Mass) ? opt.mass_by_atom.at(i) : 1.0;
      if (!(w > 0.0)) {
        throw std::runtime_error("background dielectric requires positive entity weights; check mass_field/topology masses");
      }
      wsum += w;
      rx += w * xu[i];
      ry += w * yu[i];
      rz += w * zu[i];
      qsum += q[i];
    }
    if (!(wsum > 0.0)) {
      throw std::runtime_error("background dielectric encountered an entity with non-positive total weight");
    }

    const double cx = rx / wsum;
    const double cy = ry / wsum;
    const double cz = rz / wsum;
    if (opt.charge_model == ChargeModel::NeutralOnly && std::abs(qsum) > opt.entity_charge_tolerance) {
      throw std::runtime_error(
          "background dielectric in neutral_only mode requires approximately neutral grouped entities; entity id "
          + std::to_string(entities.chain_id[c]) + " has net charge " + dstr(qsum));
    }

    std::size_t bin = 0;
    {
      double s = 0.0;
      if (slab_tf) {
        s = slab_tf->point_to_domain_fraction(frame.box, cx, cy, cz);
      } else {
        const double L = measure_ext::axis_length(frame.box, axis);
        s = measure_ext::primary_axis_coord(frame.box, cx, cy, cz, axis) / L;
      }
      s = std::clamp(s, 0.0, 1.0);
      bin = static_cast<std::size_t>(std::floor(s * static_cast<double>(n_bins)));
      if (bin >= n_bins) bin = n_bins - 1;
    }

    std::vector<double> mux(ns, 0.0), muy(ns, 0.0), muz(ns, 0.0);
    for (const std::size_t i : atoms) {
      const int sid = species.atom_species.at(i);
      if (sid < 0) continue;
      mux[sid] += q[i] * (xu[i] - cx);
      muy[sid] += q[i] * (yu[i] - cy);
      muz[sid] += q[i] * (zu[i] - cz);
    }
    for (std::size_t s = 0; s < ns; ++s) {
      out.total_mx[s] += mux[s];
      out.total_my[s] += muy[s];
      out.total_mz[s] += muz[s];
      out.bin_mx[s][bin] += mux[s];
      out.bin_my[s][bin] += muy[s];
      out.bin_mz[s][bin] += muz[s];
    }
  }

  return out;
}

class SlabBackgroundDielectricMeasure final : public IMeasure {
public:
  struct Options {
    DielectricCommonOptions common;
    Axis1D axis = Axis1D::Z;
    std::size_t n_bins = 0;
    SlabTransformOptions slab;
  };

  SlabBackgroundDielectricMeasure(std::string type_name,
                                  std::string instance_name,
                                  std::string output_path,
                                  SelectionView sel,
                                  std::optional<SelectionView> anchor_sel,
                                  SelectedChains entities,
                                  SpeciesDecomposition species,
                                  Options opt)
      : type_name_(std::move(type_name)),
        instance_name_(std::move(instance_name)),
        output_path_(std::move(output_path)),
        entities_(std::move(entities)),
        species_(std::move(species)),
        opt_(std::move(opt)) {
    sel_name_owned_ = std::string(sel.name);
    sel_idx_owned_.assign(sel.idx.begin(), sel.idx.end());
    sel_ = SelectionView{sel_name_owned_, std::span<const std::size_t>(sel_idx_owned_.data(), sel_idx_owned_.size())};
    if (anchor_sel.has_value()) {
      anchor_name_owned_ = std::string(anchor_sel->name);
      anchor_idx_owned_.assign(anchor_sel->idx.begin(), anchor_sel->idx.end());
      anchor_sel_ = SelectionView{anchor_name_owned_, std::span<const std::size_t>(anchor_idx_owned_.data(), anchor_idx_owned_.size())};
    }
    validate_common_();
    if (entities_.n_chains() == 0) throw std::runtime_error(type_name_ + ": no grouped entities are available");
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
    od.x_axis = "bin";
    od.x_unit = "fractional_slab_coordinate";
    od.columns = {"kind", "is_background", "bin", "coord_lo_frac", "coord_hi_frac", "coord_center_frac",
                  "coord_center_mean", "cov_parallel", "cov_perp", "chi_parallel_contrib", "chi_perp_contrib",
                  "frac_parallel", "frac_perp", "eps_parallel_bg", "eps_perp_bg", "n_frames"};
    md.outputs.push_back(std::move(od));

    md.params["temperature"] = dstr(opt_.common.temperature);
    md.params["charge_field"] = opt_.common.charge_field;
    md.params["length_unit"] = opt_.common.length_unit;
    md.params["charge_unit"] = opt_.common.charge_unit;
    md.params["frame_start"] = std::to_string(opt_.common.frame_start);
    md.params["frame_end"] = std::to_string(opt_.common.frame_end);
    md.params["axis"] = axis1d_name(opt_.axis);
    md.params["n_bins"] = std::to_string(opt_.n_bins);
    md.params["n_entities"] = std::to_string(entities_.n_chains());
    md.params["charge_model"] = charge_model_name(opt_.common.charge_model);
    md.params["entity_reference"] = entity_reference_name(opt_.common.entity_reference);
    md.params["assigned_atoms"] = std::to_string(species_.assigned_atoms);
    md.params["unassigned_atoms"] = std::to_string(species_.unassigned_atoms);
    md.params["allow_partial_species_cover"] = species_.allow_partial_cover ? "true" : "false";
    md.params["fraction_tolerance"] = dstr(kContributionFractionTolerance);
    md.params["slab_align_recenter"] = opt_.slab.slab_align_recenter ? "true" : "false";
    md.params["halfcell_fold"] = opt_.slab.halfcell_fold ? "true" : "false";
    md.params["target_center_frac"] = dstr(opt_.slab.target_center_frac);
    return md;
  }

  void on_start(const Frame& first_frame) override {
    (void)first_frame;
    started_ = true;
    if (!opt_.common.dry_run) flush_partial();
  }

  void on_frame(const Frame& frame, std::size_t frame_index) override {
    if (!frame_in_range(frame_index, opt_.common.frame_start, opt_.common.frame_end)) return;
    const SlabFrameTransform tf = make_slab_frame_transform(frame, anchor_sel_ ? &(*anchor_sel_) : nullptr, opt_.slab);
    frames_.push_back(compute_background_dielectric_frame(frame, entities_, opt_.common, &tf, opt_.axis,
                                                          opt_.n_bins, species_));
  }

  void flush_partial() override {
    if (opt_.common.dry_run || !started_) return;
    util::atomic_write_text(output_path_, [&](std::ostream& ofs) {
      write_header_(ofs);
      if (frames_.empty()) return;

      const double lscale = length_to_m(opt_.common.length_unit);
      const double qscale = charge_to_c(opt_.common.charge_unit);
      const double dscale = lscale * qscale;
      const double vscale = lscale * lscale * lscale;
      const double invn = 1.0 / static_cast<double>(frames_.size());
      const std::size_t ns = species_.labels.size();

      double mean_domain_length = 0.0;
      for (const auto& f : frames_) mean_domain_length += f.domain_length;
      mean_domain_length *= invn * lscale;

      for (std::size_t b = 0; b < opt_.n_bins; ++b) {
        double sum_total_bg_parallel = 0.0, sum_total_bg_perp = 0.0;
        double sum_bin_bg_parallel = 0.0, sum_bin_bg_perp = 0.0;
        double sum_bg_cross_parallel = 0.0, sum_bg_cross_perp = 0.0;
        double sum_vbin = 0.0;

        std::vector<double> sum_species_bin_parallel(ns, 0.0), sum_species_bin_perp(ns, 0.0);
        std::vector<double> sum_species_cross_parallel(ns, 0.0), sum_species_cross_perp(ns, 0.0);

        for (const auto& f : frames_) {
          std::vector<double> total_parallel(ns, 0.0), total_perp(ns, 0.0);
          std::vector<double> bin_parallel(ns, 0.0), bin_perp(ns, 0.0);

          double total_bg_parallel = 0.0, total_bg_perp = 0.0;
          double bin_bg_parallel = 0.0, bin_bg_perp = 0.0;

          for (std::size_t s = 0; s < ns; ++s) {
            const double mx = f.total_mx[s] * dscale;
            const double my = f.total_my[s] * dscale;
            const double mz = f.total_mz[s] * dscale;
            const double bmx = f.bin_mx[s][b] * dscale;
            const double bmy = f.bin_my[s][b] * dscale;
            const double bmz = f.bin_mz[s][b] * dscale;

            total_parallel[s] = parallel_average_component(opt_.axis, mx, my, mz);
            total_perp[s] = normal_component(opt_.axis, mx, my, mz);
            bin_parallel[s] = parallel_average_component(opt_.axis, bmx, bmy, bmz);
            bin_perp[s] = normal_component(opt_.axis, bmx, bmy, bmz);

            if (species_.in_background[s]) {
              total_bg_parallel += total_parallel[s];
              total_bg_perp += total_perp[s];
              bin_bg_parallel += bin_parallel[s];
              bin_bg_perp += bin_perp[s];
            }
          }

          const double vbin = (f.volume * vscale) / static_cast<double>(opt_.n_bins);
          sum_total_bg_parallel += total_bg_parallel;
          sum_total_bg_perp += total_bg_perp;
          sum_bin_bg_parallel += bin_bg_parallel;
          sum_bin_bg_perp += bin_bg_perp;
          sum_bg_cross_parallel += bin_bg_parallel * total_bg_parallel;
          sum_bg_cross_perp += bin_bg_perp * total_bg_perp;
          sum_vbin += vbin;

          for (std::size_t s = 0; s < ns; ++s) {
            sum_species_bin_parallel[s] += bin_parallel[s];
            sum_species_bin_perp[s] += bin_perp[s];
            sum_species_cross_parallel[s] += bin_parallel[s] * total_bg_parallel;
            sum_species_cross_perp[s] += bin_perp[s] * total_bg_perp;
          }
        }

        const double mean_total_bg_parallel = sum_total_bg_parallel * invn;
        const double mean_total_bg_perp = sum_total_bg_perp * invn;
        const double mean_bin_bg_parallel = sum_bin_bg_parallel * invn;
        const double mean_bin_bg_perp = sum_bin_bg_perp * invn;
        const double mean_vbin = sum_vbin * invn;

        const double cov_bg_parallel = sum_bg_cross_parallel * invn - mean_bin_bg_parallel * mean_total_bg_parallel;
        const double cov_bg_perp = sum_bg_cross_perp * invn - mean_bin_bg_perp * mean_total_bg_perp;

        const double chi_bg_parallel = cov_bg_parallel / (kEps0 * mean_vbin * kBoltzmann * opt_.common.temperature);
        const double chi_bg_perp = cov_bg_perp / (kEps0 * mean_vbin * kBoltzmann * opt_.common.temperature);
        const double eps_bg_parallel = 1.0 + chi_bg_parallel;
        const double denom = 1.0 - chi_bg_perp;
        const double eps_bg_perp = (std::abs(denom) > 1e-12)
                                     ? (1.0 / denom)
                                     : ((denom >= 0.0) ? 1e12 : -1e12);

        const double lo_frac = static_cast<double>(b) / static_cast<double>(opt_.n_bins);
        const double hi_frac = static_cast<double>(b + 1) / static_cast<double>(opt_.n_bins);
        const double center_frac = 0.5 * (lo_frac + hi_frac);
        const double center = center_frac * mean_domain_length;

        auto write_row = [&](const std::string& label,
                             int is_bg,
                             double cov_parallel,
                             double cov_perp,
                             double chi_parallel,
                             double chi_perp,
                             double frac_parallel,
                             double frac_perp) {
          ofs << label << ' '
              << is_bg << ' '
              << b << ' '
              << std::setprecision(17) << lo_frac << ' '
              << std::setprecision(17) << hi_frac << ' '
              << std::setprecision(17) << center_frac << ' '
              << std::setprecision(17) << center << ' '
              << std::setprecision(17) << cov_parallel << ' '
              << std::setprecision(17) << cov_perp << ' '
              << std::setprecision(17) << chi_parallel << ' '
              << std::setprecision(17) << chi_perp << ' '
              << std::setprecision(17) << frac_parallel << ' '
              << std::setprecision(17) << frac_perp << ' '
              << std::setprecision(17) << eps_bg_parallel << ' '
              << std::setprecision(17) << eps_bg_perp << ' '
              << frames_.size() << '\n';
        };

        write_row("TOTAL_BG", 1, cov_bg_parallel, cov_bg_perp, chi_bg_parallel, chi_bg_perp, 1.0, 1.0);

        for (std::size_t s = 0; s < ns; ++s) {
          const double mean_species_bin_parallel = sum_species_bin_parallel[s] * invn;
          const double mean_species_bin_perp = sum_species_bin_perp[s] * invn;
          const double cov_parallel = sum_species_cross_parallel[s] * invn
                                    - mean_species_bin_parallel * mean_total_bg_parallel;
          const double cov_perp = sum_species_cross_perp[s] * invn
                                - mean_species_bin_perp * mean_total_bg_perp;
          const double chi_parallel = cov_parallel / (kEps0 * mean_vbin * kBoltzmann * opt_.common.temperature);
          const double chi_perp = cov_perp / (kEps0 * mean_vbin * kBoltzmann * opt_.common.temperature);
          const double frac_parallel = contribution_fraction_or_nan(chi_parallel, chi_bg_parallel);
          const double frac_perp = contribution_fraction_or_nan(chi_perp, chi_bg_perp);
          write_row(species_.labels[s], species_.in_background[s] ? 1 : 0,
                    cov_parallel, cov_perp, chi_parallel, chi_perp, frac_parallel, frac_perp);
        }
      }
    });
  }

  void finalize() override {
    if (opt_.common.dry_run) return;
    flush_partial();
  }

private:
  std::string type_name_;
  std::string instance_name_;
  std::string output_path_;
  std::string sel_name_owned_;
  std::vector<std::size_t> sel_idx_owned_;
  SelectionView sel_{{}, {}};
  std::string anchor_name_owned_;
  std::vector<std::size_t> anchor_idx_owned_;
  std::optional<SelectionView> anchor_sel_;
  SelectedChains entities_;
  SpeciesDecomposition species_;
  Options opt_;
  std::vector<BackgroundDielectricFrame> frames_;
  bool started_ = false;

  void validate_common_() const {
    if (opt_.common.frame_start < 0) throw std::runtime_error(type_name_ + ": frame_start must be >= 0");
    if (opt_.common.frame_end >= 0 && opt_.common.frame_end < opt_.common.frame_start) {
      throw std::runtime_error(type_name_ + ": frame_end must be -1 or >= frame_start");
    }
    if (!(opt_.common.temperature > 0.0)) throw std::runtime_error(type_name_ + ": temperature must be > 0");
    if (opt_.n_bins == 0) throw std::runtime_error(type_name_ + ": n_bins must be >= 1");
    if (sel_.idx.empty()) throw std::runtime_error(type_name_ + ": selection is empty");
    if (species_.labels.empty()) throw std::runtime_error(type_name_ + ": no species groups were defined");
    if (opt_.common.entity_reference == EntityReference::Mass && opt_.common.mass_by_atom.empty()) {
      throw std::runtime_error(type_name_ + ": entity_reference=mass requires mass_field or topology masses");
    }
  }

  void write_header_(std::ostream& ofs) const {
    ofs << "# PILOTS: slab background dielectric profile from species-resolved local polarization-density fluctuations\n";
    ofs << "# selection: " << sel_.name << " (n=" << sel_.idx.size() << ")\n";
    ofs << "# n_entities: " << entities_.n_chains() << "\n";
    ofs << "# axis: " << axis1d_name(opt_.axis) << ", n_bins: " << opt_.n_bins << "\n";
    ofs << "# charge_field: " << opt_.common.charge_field << "\n";
    ofs << "# charge_unit: " << opt_.common.charge_unit << ", length_unit: " << opt_.common.length_unit << "\n";
    ofs << "# temperature: " << std::setprecision(17) << opt_.common.temperature << " K\n";
    ofs << "# frame_range: [" << opt_.common.frame_start << ", " << opt_.common.frame_end << "]\n";
    ofs << "# charge_model: " << charge_model_name(opt_.common.charge_model) << "\n";
    ofs << "# entity_reference: " << entity_reference_name(opt_.common.entity_reference) << "\n";
    ofs << "# slab_align_recenter: " << (opt_.slab.slab_align_recenter ? "true" : "false")
        << ", halfcell_fold: " << (opt_.slab.halfcell_fold ? "true" : "false")
        << ", target_center_frac: " << std::setprecision(17) << opt_.slab.target_center_frac << "\n";
    if (anchor_sel_) {
      ofs << "# anchor_selection: " << anchor_sel_->name << " (n=" << anchor_sel_->idx.size() << ")\n";
    }
    ofs << "# species_groups:";
    for (std::size_t s = 0; s < species_.labels.size(); ++s) {
      ofs << ' ' << species_.labels[s] << "(n=" << species_.selected_counts[s]
          << ",bg=" << (species_.in_background[s] ? "1" : "0") << ")";
    }
    ofs << "\n";
    ofs << "# assigned_atoms: " << species_.assigned_atoms
        << ", unassigned_atoms: " << species_.unassigned_atoms << "\n";
    ofs << "# allow_partial_species_cover: "
        << (species_.allow_partial_cover ? "true" : "false") << "\n";
    if (species_.unassigned_atoms > 0) {
      ofs << "# warning: species decomposition is partial; atoms outside species_groups still contribute to entity reference points and entity net-charge checks.\n";
    }
    ofs << "# fraction_tolerance: frac_parallel/frac_perp are NaN when |TOTAL_BG susceptibility| <= "
        << std::setprecision(17) << kContributionFractionTolerance << "\n";
    ofs << "# meaning: species-resolved background dielectric profile. TOTAL_BG is the background-polarization dielectric built from the selected background_groups. Species rows report additive susceptibility contributions (or cross-correlated non-background contributions) with respect to TOTAL_BG.\n";
    ofs << "# columns: kind  is_background  bin  coord_lo_frac  coord_hi_frac  coord_center_frac  coord_center_mean  cov_parallel  cov_perp  chi_parallel_contrib  chi_perp_contrib  frac_parallel  frac_perp  eps_parallel_bg  eps_perp_bg  n_frames\n";
  }
};

void append_slab_background_dielectric_caps(const IniConfig& cfg,
                                            const std::string& section,
                                            const MeasureBuildEnv& env,
                                            MeasureCapabilities& caps) {
  append_dielectric_caps(cfg, section, env, caps);
  if (!cfg.has_key(section, "species_groups")) return;
  const auto groups = get_required_name_list(cfg, section, "species_groups");
  for (const auto& g : groups) caps.group_refs.push_back(g);
  if (cfg.has_key(section, "background_groups")) {
    const auto bg = get_required_name_list(cfg, section, "background_groups");
    for (const auto& g : bg) caps.group_refs.push_back(g);
  }
}

MeasureCapabilities slab_background_dielectric_caps(const IniConfig& cfg,
                                                    const std::string& section,
                                                    const std::string& instance,
                                                    const MeasureBuildEnv& env) {
  (void)instance;
  MeasureCapabilities caps;
  append_slab_background_dielectric_caps(cfg, section, env, caps);
  return caps;
}

std::unique_ptr<IMeasure> slab_background_dielectric_create_impl(const std::string& type_name,
                                                                 const IniConfig& cfg,
                                                                 const std::string& section,
                                                                 const std::string& instance,
                                                                 const MeasureBuildEnv& env,
                                                                 const SystemContext& sysctx) {
  if (!env.first_frame || !env.selection_provider) {
    throw std::runtime_error(type_name + " factory: missing first_frame or SelectionProvider");
  }
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, type_name);
  auto anchor = build_optional_anchor_selection(cfg, section, *env.selection_provider, frame0, sel,
                                                slab_transform_requested(cfg, section), type_name,
                                                std::optional<std::string>(group),
                                                std::optional<std::string>(topo),
                                                std::optional<std::string>(comb));
  const auto entity_ids = entity_id_per_atom_from_config(cfg, section, frame0, sysctx);
  SelectedChains entities = build_selected_chains(sel, entity_ids);
  SpeciesDecomposition species = build_species_decomposition(cfg, section, frame0, *env.selection_provider, sel, type_name);

  const fs::path out = resolve_measure_output_path(cfg, section, env, "output", type_name + ".dat");
  SlabBackgroundDielectricMeasure::Options opt;
  opt.common = common_options_from_config(cfg, section, frame0, sysctx, ChargeModel::FragmentInternal, env.dry_run);
  opt.axis = parse_axis1d(cfg.get_string(section, "axis", std::optional<std::string>("z")));
  opt.n_bins = static_cast<std::size_t>(cfg.get_int64(section, "n_bins", std::optional<std::int64_t>(100)));
  opt.slab.axis = opt.axis;
  opt.slab.slab_align_recenter = cfg.get_bool(section, "slab_align_recenter", std::optional<bool>(false));
  opt.slab.halfcell_fold = cfg.get_bool(section, "halfcell_fold", std::optional<bool>(false));
  opt.slab.target_center_frac = cfg.get_double(section, "target_center_frac", std::optional<double>(0.5));
  return std::make_unique<SlabBackgroundDielectricMeasure>(type_name, instance, out.string(), sel, anchor,
                                                           std::move(entities), std::move(species), std::move(opt));
}

std::unique_ptr<IMeasure> slab_bg_dielectric_create(const IniConfig& cfg,
                                                    const std::string& section,
                                                    const std::string& instance,
                                                    const MeasureBuildEnv& env,
                                                    const SystemContext& sysctx) {
  return slab_background_dielectric_create_impl("slab_bg_dielectric", cfg, section, instance, env, sysctx);
}

void append_dielectric_caps(const IniConfig& cfg,
                            const std::string& section,
                            const MeasureBuildEnv& env,
                            MeasureCapabilities& caps) {
  caps.selection_policy = SelectionPolicy::RequireStatic;
  caps.requires_identity_consistent = true;
  caps.requires_dfields = {"xu", "yu", "zu", cfg.get_string(section, "charge_field", std::optional<std::string>("q"))};
  caps.scale = ScaleCompatibility{true, true, true};

  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  (void)topo;
  caps.group_refs.push_back(group);

  const bool want_anchor = slab_transform_requested(cfg, section);
  if (want_anchor) {
    caps.group_refs.push_back(cfg.get_string(section, "anchor_group", std::optional<std::string>(group)));
  }

  if (cfg.has_key(section, "entity_id_field")) {
    append_integer_like_field_cap(caps, cfg.get_string(section, "entity_id_field"));
  } else if (env.first_frame && env.first_frame->has_mol) {
    caps.requires_i64fields.push_back("mol");
  } else {
    caps.requires_topology_sections.push_back("bonds");
  }

  const EntityReference ref = parse_entity_reference(cfg.get_string(section, "entity_reference", std::optional<std::string>("geometric")));
  if (ref == EntityReference::Mass) {
    if (cfg.has_key(section, "mass_field")) {
      caps.requires_dfields.push_back(cfg.get_string(section, "mass_field"));
    } else {
      caps.requires_intfields.push_back("type");
      caps.requires_topology_sections.push_back("masses");
    }
  }
}

DielectricCommonOptions common_options_from_config(const IniConfig& cfg,
                                                    const std::string& section,
                                                    const Frame& frame0,
                                                    const SystemContext& sysctx,
                                                    ChargeModel default_model,
                                                    bool dry_run) {
  DielectricCommonOptions opt;
  opt.frame_start = cfg.get_int64(section, "frame_start", std::optional<std::int64_t>(0));
  opt.frame_end = cfg.get_int64(section, "frame_end", std::optional<std::int64_t>(-1));
  opt.charge_field = cfg.get_string(section, "charge_field", std::optional<std::string>("q"));
  opt.length_unit = cfg.get_string(section, "length_unit", std::optional<std::string>("angstrom"));
  opt.charge_unit = cfg.get_string(section, "charge_unit", std::optional<std::string>("e"));
  opt.temperature = cfg.get_double(section, "temperature");
  opt.entity_charge_tolerance = cfg.get_double(section, "entity_charge_tolerance", std::optional<double>(1e-6));
  opt.charge_model = charge_model_from_config(cfg, section, default_model);
  opt.entity_reference = parse_entity_reference(cfg.get_string(section, "entity_reference", std::optional<std::string>("geometric")));
  opt.dry_run = dry_run;
  if (opt.entity_reference == EntityReference::Mass) {
    opt.mass_by_atom = mass_by_atom_from_config(cfg, section, frame0, sysctx);
    if (opt.mass_by_atom.size() != frame0.natoms) {
      throw std::runtime_error("dielectric measure: mass_by_atom size mismatch");
    }
  }
  return opt;
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

MeasureCapabilities slab_dielectric_caps(const IniConfig& cfg,
                                         const std::string& section,
                                         const std::string& instance,
                                         const MeasureBuildEnv& env) {
  (void)instance;
  MeasureCapabilities caps;
  append_dielectric_caps(cfg, section, env, caps);
  return caps;
}

std::unique_ptr<IMeasure> bulk_dielectric_create_impl(const std::string& type_name,
                                                      ChargeModel default_model,
                                                      const IniConfig& cfg,
                                                      const std::string& section,
                                                      const std::string& instance,
                                                      const MeasureBuildEnv& env,
                                                      const SystemContext& sysctx) {
  if (!env.first_frame || !env.selection_provider) {
    throw std::runtime_error(type_name + " factory: missing first_frame or SelectionProvider");
  }
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, type_name);
  const auto entity_ids = entity_id_per_atom_from_config(cfg, section, frame0, sysctx);
  SelectedChains entities = build_selected_chains(sel, entity_ids);

  const fs::path out = resolve_measure_output_path(cfg, section, env, "output", type_name + ".dat");
  DielectricCommonOptions common = common_options_from_config(cfg, section, frame0, sysctx, default_model, env.dry_run);
  return std::make_unique<BulkDielectricMeasure>(type_name, instance, out.string(), sel, std::move(entities), std::move(common));
}

std::unique_ptr<IMeasure> slab_dielectric_create_impl(const std::string& type_name,
                                                      ChargeModel default_model,
                                                      const IniConfig& cfg,
                                                      const std::string& section,
                                                      const std::string& instance,
                                                      const MeasureBuildEnv& env,
                                                      const SystemContext& sysctx) {
  if (!env.first_frame || !env.selection_provider) {
    throw std::runtime_error(type_name + " factory: missing first_frame or SelectionProvider");
  }
  const Frame& frame0 = *env.first_frame;
  const std::string group = cfg.get_string(section, "group", std::optional<std::string>("all"));
  const std::string topo = cfg.get_string(section, "topo_group", std::optional<std::string>("all"));
  const std::string comb = cfg.get_string(section, "combine", std::optional<std::string>("A&T"));
  SelectionView sel = get_static_combined_view(*env.selection_provider, frame0, group, topo, comb, type_name);
  auto anchor = build_optional_anchor_selection(cfg, section, *env.selection_provider, frame0, sel,
                                                slab_transform_requested(cfg, section), type_name,
                                                std::optional<std::string>(group),
                                                std::optional<std::string>(topo),
                                                std::optional<std::string>(comb));
  const auto entity_ids = entity_id_per_atom_from_config(cfg, section, frame0, sysctx);
  SelectedChains entities = build_selected_chains(sel, entity_ids);

  const fs::path out = resolve_measure_output_path(cfg, section, env, "output", type_name + ".dat");

  SlabDielectricMeasure::Options opt;
  opt.common = common_options_from_config(cfg, section, frame0, sysctx, default_model, env.dry_run);
  opt.axis = parse_axis1d(cfg.get_string(section, "axis", std::optional<std::string>("z")));
  opt.n_bins = static_cast<std::size_t>(cfg.get_int64(section, "n_bins", std::optional<std::int64_t>(100)));
  opt.slab.axis = opt.axis;
  opt.slab.slab_align_recenter = cfg.get_bool(section, "slab_align_recenter", std::optional<bool>(false));
  opt.slab.halfcell_fold = cfg.get_bool(section, "halfcell_fold", std::optional<bool>(false));
  opt.slab.target_center_frac = cfg.get_double(section, "target_center_frac", std::optional<double>(0.5));

  return std::make_unique<SlabDielectricMeasure>(type_name, instance, out.string(), sel, anchor, std::move(entities), std::move(opt));
}

std::unique_ptr<IMeasure> bulk_dielectric_create(const IniConfig& cfg,
                                                 const std::string& section,
                                                 const std::string& instance,
                                                 const MeasureBuildEnv& env,
                                                 const SystemContext& sysctx) {
  return bulk_dielectric_create_impl("bulk_dielectric", ChargeModel::NeutralOnly, cfg, section, instance, env, sysctx);
}

std::unique_ptr<IMeasure> slab_dielectric_create(const IniConfig& cfg,
                                                 const std::string& section,
                                                 const std::string& instance,
                                                 const MeasureBuildEnv& env,
                                                 const SystemContext& sysctx) {
  return slab_dielectric_create_impl("slab_dielectric", ChargeModel::NeutralOnly, cfg, section, instance, env, sysctx);
}

std::unique_ptr<IMeasure> bulk_fragment_dielectric_create(const IniConfig& cfg,
                                                          const std::string& section,
                                                          const std::string& instance,
                                                          const MeasureBuildEnv& env,
                                                          const SystemContext& sysctx) {
  return bulk_dielectric_create_impl("bulk_fragment_dielectric", ChargeModel::FragmentInternal, cfg, section, instance, env, sysctx);
}

std::unique_ptr<IMeasure> slab_fragment_dielectric_create(const IniConfig& cfg,
                                                          const std::string& section,
                                                          const std::string& instance,
                                                          const MeasureBuildEnv& env,
                                                          const SystemContext& sysctx) {
  return slab_dielectric_create_impl("slab_fragment_dielectric", ChargeModel::FragmentInternal, cfg, section, instance, env, sysctx);
}

static MeasureRegistrar g_register_bulk_dielectric("bulk_dielectric", &bulk_dielectric_caps, &bulk_dielectric_create);
static MeasureRegistrar g_register_slab_dielectric("slab_dielectric", &slab_dielectric_caps, &slab_dielectric_create);
static MeasureRegistrar g_register_slab_dielectric_alias("layered_dielectric", &slab_dielectric_caps, &slab_dielectric_create);
static MeasureRegistrar g_register_bulk_fragment_dielectric("bulk_fragment_dielectric", &bulk_dielectric_caps, &bulk_fragment_dielectric_create);
static MeasureRegistrar g_register_slab_fragment_dielectric("slab_fragment_dielectric", &slab_dielectric_caps, &slab_fragment_dielectric_create);
static MeasureRegistrar g_register_layered_fragment_dielectric_alias("layered_fragment_dielectric", &slab_dielectric_caps, &slab_fragment_dielectric_create);
static MeasureRegistrar g_register_slab_bg_dielectric("slab_bg_dielectric", &slab_background_dielectric_caps, &slab_bg_dielectric_create);
static MeasureRegistrar g_register_layered_bg_dielectric_alias("layered_bg_dielectric", &slab_background_dielectric_caps, &slab_bg_dielectric_create);

} // namespace
} // namespace pilots
