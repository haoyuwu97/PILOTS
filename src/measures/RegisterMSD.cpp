#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

#include "MeasureCommon.hpp"
#include "pilots/correlate/CorrelatorFactory.hpp"
#include "pilots/measures/MSDMeasure.hpp"
#include "pilots/measures/MeasureRegistry.hpp"

namespace fs = std::filesystem;

namespace pilots {
namespace {

using measure_ext::get_static_combined_view;
using measure_ext::get_static_group_view;
using measure_ext::parse_diag_mask;
using measure_ext::resolve_exact_frame_end;
using measure_ext::resolve_path;

MeasureCapabilities msd_caps(const IniConfig& cfg,
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

std::unique_ptr<IMeasure> msd_create(const IniConfig& cfg,
                                    const std::string& section,
                                    const std::string& instance,
                                    const MeasureBuildEnv& env,
                                    const SystemContext& sysctx) {
  (void)sysctx;

  if (!env.first_frame) throw std::runtime_error("MSD factory: first_frame is null");
  if (!env.selection_provider) {
    throw std::runtime_error("MSD factory: SelectionProvider is missing (Runner must build it)");
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
  if (frame_end >= 0 && frame_end < frame_start) throw std::runtime_error("frame_end must be -1 or >= frame_start");

  const fs::path out_path = measure_ext::resolve_measure_output_path(cfg, section, env, "output", "msd.dat");

  const std::int64_t frame_end_eff = resolve_exact_frame_end(
      corr_spec, env.follow, frame_start, frame_end, env.input_path);

  const int diag_mask = parse_diag_mask(comp_s);

  // Selections (AtomGroup + TopoGroup + combine).
  // MSD is scientifically meaningful only for identity-consistent, static selections.
  SelectionView sel = get_static_combined_view(
      *env.selection_provider, frame0, group_ref, topo_group_ref, combine_expr, "MSD");

  SelectionView drift_sel;
  if (remove_drift) {
    drift_sel = get_static_group_view(
        *env.selection_provider, frame0, drift_group_ref, "MSD");
  } else {
    drift_sel = get_static_group_view(*env.selection_provider, frame0, "all", "MSD");
  }

  MSDMeasure::Options opt;
  opt.frame_start = frame_start;
  opt.frame_end = frame_end_eff;
  opt.diag_mask = diag_mask;
  opt.remove_drift = remove_drift;
  opt.corr = corr_spec;
  opt.dry_run = env.dry_run;

  return std::make_unique<MSDMeasure>(instance, out_path.string(), sel, drift_sel, opt);
}

static MeasureRegistrar g_register_msd("msd", &msd_caps, &msd_create);

} // namespace
} // namespace pilots
