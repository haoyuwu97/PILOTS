#include "pilots/app/Runner.hpp"

#include <cctype>
#include <chrono>
#include <algorithm>
#include <csignal>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#if PILOTS_HAS_OPENMP
#include <omp.h>
#endif

#include "pilots/core/MoleculeIndex.hpp"
#include "pilots/core/SystemContext.hpp"
#include "pilots/io/LammpsDumpReader.hpp"
#include "pilots/io/TopologyReaders.hpp"
#include "pilots/measures/IMeasure.hpp"
#include "pilots/measures/MeasureRegistry.hpp"
#include "pilots/output/ResultsIndex.hpp"
#include "pilots/select/GroupRegistry.hpp"
#include "pilots/select/TopoGroupRegistry.hpp"
#include "pilots/select/SelectionProvider.hpp"
#include "pilots/topology/Topology.hpp"
#include "pilots/topology/BondGraph.hpp"
#include "pilots/alg/mapping/BeadGraphBuilder.hpp"
#include "pilots/alg/mapping/BeadMapping.hpp"
#include "pilots/alg/mapping/MappingSpec.hpp"
#include "pilots/alg/polymer/PolymerClassifier.hpp"
#include "pilots/alg/index/ChainIndex.hpp"
#include "pilots/util/AtomicFile.hpp"
#include "pilots/util/BinaryIO.hpp"
#include "pilots/util/Hash.hpp"
#include "pilots/util/Timer.hpp"

namespace fs = std::filesystem;

namespace {
static volatile std::sig_atomic_t g_stop = 0;
static void handle_signal(int) { g_stop = 1; }

fs::path resolve_path(const fs::path& base_dir, const std::string& p) {
  fs::path path(p);
  if (path.is_absolute()) return path;
  return (base_dir / path).lexically_normal();
}

struct ResumeMeta {
  std::uint64_t config_hash = 0;
  std::string input_path;
  // Hash of mapping spec (K6), to ensure resume uses identical bead definition.
  std::string mapping_spec_hash_hex;
  std::uint64_t reader_offset = 0;
  std::uint64_t frame_idx = 0; // next frame index to process
  double wall_seconds = 0.0;
  double reader_seconds = 0.0;
  double group_setup_seconds = 0.0;
  // v2+ topology meta (auditable; used for sanity checks on resume)
  std::string topology_file;
  std::string topology_format;
  std::string topology_hash_fnv1a64_hex;
  std::vector<std::string> topology_loaded_sections;
  bool topology_derive_mol_from_bonds = false;
  std::vector<pilots::output::MeasureProfiling> mp;
};

inline std::string hex_u64(std::uint64_t h) {
  std::ostringstream oss;
  oss << std::hex << std::setw(16) << std::setfill('0') << h;
  return oss.str();
}

inline std::int64_t file_time_to_epoch_seconds(fs::file_time_type t) {
  // Best-effort conversion from filesystem clock to system_clock.
  // This avoids relying on C++20 chrono extensions while keeping results.json auditable.
  using namespace std::chrono;
  const auto now_fs = fs::file_time_type::clock::now();
  const auto now_sys = system_clock::now();
  const auto sys_time = time_point_cast<system_clock::duration>(t - now_fs + now_sys);
  return duration_cast<seconds>(sys_time.time_since_epoch()).count();
}

inline pilots::output::FileFingerprint make_fingerprint(const fs::path& p,
                                                        bool compute_hash,
                                                        std::optional<std::uint64_t> known_hash = std::nullopt) {
  pilots::output::FileFingerprint fp;
  fp.path = p.string();
  try {
    fp.size_bytes = fs::exists(p) ? static_cast<std::uint64_t>(fs::file_size(p)) : 0ull;
  } catch (...) {
    fp.size_bytes = 0ull;
  }
  try {
    fp.mtime_epoch_s = fs::exists(p) ? file_time_to_epoch_seconds(fs::last_write_time(p)) : 0;
  } catch (...) {
    fp.mtime_epoch_s = 0;
  }

  if (compute_hash) {
    fp.hash_kind = "fnv1a64";
    const std::uint64_t h = known_hash ? *known_hash : pilots::fnv1a64_file(p.string());
    fp.hash_fnv1a64_hex = hex_u64(h);
    fp.hash_computed = true;
  } else {
    fp.hash_kind = "none";
    fp.hash_fnv1a64_hex.clear();
    fp.hash_computed = false;
  }
  return fp;
}

ResumeMeta load_checkpoint_into(const fs::path& ckpt_path,
                               std::uint64_t expected_config_hash,
                               const fs::path& expected_input_path,
                               const std::string& expected_mapping_spec_hash_hex,
                               const std::vector<std::unique_ptr<pilots::IMeasure>>& measures,
                               std::vector<pilots::output::MeasureProfiling>& mp) {
  std::ifstream ifs(ckpt_path, std::ios::binary);
  if (!ifs) throw std::runtime_error("failed to open checkpoint: " + ckpt_path.string());

  pilots::util::require_magic(ifs, "PILOTSCHKPT");
  pilots::util::BinaryReader r(ifs);
  const std::uint32_t ver = r.read_u32();
  if (ver != 1 && ver != 2 && ver != 3) throw std::runtime_error("checkpoint version not supported");

  (void)r.read_string(); // pilots_version at save (informational)

  ResumeMeta meta;
  meta.config_hash = r.read_u64();
  if (meta.config_hash != expected_config_hash) {
    throw std::runtime_error("checkpoint config_hash mismatch (wrong resume_from or config changed)");
  }

  meta.input_path = r.read_string();
  if (fs::path(meta.input_path) != expected_input_path) {
    throw std::runtime_error("checkpoint input_path mismatch (wrong resume_from or config/input changed)\n"
                             "  checkpoint: " + meta.input_path + "\n"
                             "  current:    " + expected_input_path.string());
  }

  // Mapping spec sanity (K6): we require mapping spec hash to match when mapping is enabled.
  // For older checkpoints (ver<3), we refuse to resume into a mapping-enabled run because
  // bead definitions cannot be validated.
  if (ver >= 3) {
    meta.mapping_spec_hash_hex = r.read_string();
  } else {
    meta.mapping_spec_hash_hex.clear();
  }

  const bool want_mapping = !expected_mapping_spec_hash_hex.empty();
  const bool have_mapping_meta = !meta.mapping_spec_hash_hex.empty();
  if (want_mapping) {
    if (!have_mapping_meta) {
      throw std::runtime_error("checkpoint lacks mapping spec metadata (ver<3); cannot safely resume with mapping enabled");
    }
    if (meta.mapping_spec_hash_hex != expected_mapping_spec_hash_hex) {
      throw std::runtime_error("checkpoint mapping spec hash mismatch (mapping config or file changed)");
    }
  } else {
    // If checkpoint was created with mapping, require current run to also specify mapping.
    if (have_mapping_meta) {
      throw std::runtime_error("checkpoint was created with mapping enabled, but current run has mapping disabled");
    }
  }

  meta.reader_offset = r.read_u64();
  meta.frame_idx = r.read_u64();
  meta.wall_seconds = r.read_f64();
  meta.reader_seconds = r.read_f64();
  meta.group_setup_seconds = r.read_f64();

  if (ver >= 2) {
    meta.topology_file = r.read_string();
    meta.topology_format = r.read_string();
    meta.topology_hash_fnv1a64_hex = r.read_string();
    const std::uint32_t nsec = r.read_u32();
    meta.topology_loaded_sections.clear();
    meta.topology_loaded_sections.reserve(nsec);
    for (std::uint32_t i = 0; i < nsec; ++i) {
      meta.topology_loaded_sections.push_back(r.read_string());
    }
    meta.topology_derive_mol_from_bonds = (r.read_u8() != 0);
  }

  const std::uint64_t n_meas = r.read_u64();
  if (n_meas != static_cast<std::uint64_t>(measures.size())) {
    throw std::runtime_error("checkpoint measure count mismatch");
  }

  meta.mp.resize(measures.size());
  for (std::size_t i = 0; i < measures.size(); ++i) {
    const std::string inst = r.read_string();
    if (inst != measures[i]->instance_name()) {
      throw std::runtime_error("checkpoint measure order/name mismatch: expected '" + measures[i]->instance_name() + "' got '" + inst + "'");
    }
    meta.mp[i].on_start_s = r.read_f64();
    meta.mp[i].on_frame_s = r.read_f64();
    meta.mp[i].finalize_s = r.read_f64();
    meta.mp[i].frames = static_cast<std::size_t>(r.read_u64());

    measures[i]->load_state(ifs);
  }

  pilots::util::require_magic(ifs, "PILOTSEND");

  // merge profiling into current mp vector
  if (mp.size() != measures.size()) mp.assign(measures.size(), pilots::output::MeasureProfiling{});
  for (std::size_t i = 0; i < measures.size(); ++i) {
    mp[i].on_start_s += meta.mp[i].on_start_s;
    mp[i].on_frame_s += meta.mp[i].on_frame_s;
    mp[i].finalize_s += meta.mp[i].finalize_s;
    mp[i].frames += meta.mp[i].frames;
  }

  return meta;
}

void save_checkpoint_from(const fs::path& ckpt_path,
                          std::uint64_t config_hash,
                          const fs::path& input_path,
                          const std::string& mapping_spec_hash_hex,
                          std::uint64_t reader_offset,
                          std::uint64_t frame_idx,
                          double wall_seconds,
                          double reader_seconds,
                          double group_setup_seconds,
                          const std::string& topology_file,
                          const std::string& topology_format,
                          const std::string& topology_hash_fnv1a64_hex,
                          const std::vector<std::string>& topology_loaded_sections,
                          bool topology_derive_mol_from_bonds,
                          const std::vector<std::unique_ptr<pilots::IMeasure>>& measures,
                          const std::vector<pilots::output::MeasureProfiling>& mp) {
  pilots::util::atomic_write_binary(ckpt_path, [&](std::ostream& os) {
    pilots::util::write_magic(os, "PILOTSCHKPT");
    pilots::util::BinaryWriter w(os);
    w.write_u32(3);
    w.write_string(PILOTS_VERSION_STR);

    w.write_u64(config_hash);
    w.write_string(input_path.string());
    w.write_string(mapping_spec_hash_hex);
    w.write_u64(reader_offset);
    w.write_u64(frame_idx);
    w.write_f64(wall_seconds);
    w.write_f64(reader_seconds);
    w.write_f64(group_setup_seconds);

    // v2 meta
    w.write_string(topology_file);
    w.write_string(topology_format);
    w.write_string(topology_hash_fnv1a64_hex);
    w.write_u32(static_cast<std::uint32_t>(topology_loaded_sections.size()));
    for (const auto& s : topology_loaded_sections) {
      w.write_string(s);
    }
    w.write_u8(static_cast<std::uint8_t>(topology_derive_mol_from_bonds ? 1 : 0));

    w.write_u64(static_cast<std::uint64_t>(measures.size()));
    for (std::size_t i = 0; i < measures.size(); ++i) {
      w.write_string(measures[i]->instance_name());
      w.write_f64(mp[i].on_start_s);
      w.write_f64(mp[i].on_frame_s);
      w.write_f64(mp[i].finalize_s);
      w.write_u64(static_cast<std::uint64_t>(mp[i].frames));

      measures[i]->save_state(os);
    }

    pilots::util::write_magic(os, "PILOTSEND");
  });
}

std::vector<std::string> parse_csv_list(std::string s) {
  std::vector<std::string> out;
  std::string cur;
  auto flush = [&]() {
    // trim
    std::size_t b = 0;
    while (b < cur.size() && (cur[b] == ' ' || cur[b] == '\t' || cur[b] == '\r' || cur[b] == '\n')) ++b;
    std::size_t e = cur.size();
    while (e > b && (cur[e - 1] == ' ' || cur[e - 1] == '\t' || cur[e - 1] == '\r' || cur[e - 1] == '\n')) --e;
    if (e > b) out.push_back(cur.substr(b, e - b));
    cur.clear();
  };

  for (char ch : s) {
    if (ch == ',') {
      flush();
    } else {
      cur.push_back(ch);
    }
  }
  flush();

  // normalize lower-case for section names
  for (auto& x : out) {
    for (auto& c : x) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
  }
  // drop empties
  out.erase(std::remove_if(out.begin(), out.end(), [](const std::string& x) { return x.empty(); }), out.end());
  // unique
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

using TopologyLoadMask = pilots::TopologySectionMask;

TopologyLoadMask to_topology_mask(const std::vector<std::string>& sections) {
  TopologyLoadMask m;
  for (const auto& s : sections) {
    if (s == "bonds" || s == "bond") m.bonds = true;
    else if (s == "masses" || s == "mass") m.masses = true;
    else if (s == "angles" || s == "angle") m.angles = true;
    else if (s == "dihedrals" || s == "dihedral") m.dihedrals = true;
    else if (s == "impropers" || s == "improper") m.impropers = true;
    else {
      throw std::runtime_error("unknown topology section name: '" + s + "'");
    }
  }
  return m;
}

} // namespace

namespace pilots {

Runner::Runner(const IniConfig& cfg) : cfg_(cfg) {}

int Runner::run() {
  return run_impl_(false);
}

int Runner::validate_config() {
  return run_impl_(true);
}

int Runner::run_impl_(bool validate_only) {
  WallTimer total_timer;

  // --- General ---
  const fs::path cfg_path = cfg_.file_path();
  const fs::path cfg_dir = cfg_.base_dir();

  const std::string input_s = cfg_.get_string("general", "input");
  const double dt = cfg_.get_double("general", "dt");
  if (!(dt > 0.0)) throw std::runtime_error("dt must be > 0");

  const std::string outdir_s = cfg_.get_string("general", "output_dir", std::optional<std::string>("./out"));
  const bool print_profile = cfg_.get_bool("general", "profile", std::optional<bool>(true));
  const std::string results_json_name = cfg_.get_string("general", "results_json", std::optional<std::string>("results.json"));

  const bool follow = cfg_.get_bool("general", "follow", std::optional<bool>(false));
  const int follow_poll_ms = static_cast<int>(cfg_.get_int64("general", "follow_poll_ms", std::optional<std::int64_t>(200)));
  const bool allow_header_change = cfg_.get_bool("general", "allow_header_change", std::optional<bool>(false));

  const std::size_t flush_every_frames = cfg_.get_size("general", "flush_every_frames", std::optional<std::size_t>(0));
  const double flush_every_seconds = cfg_.get_double("general", "flush_every_seconds", std::optional<double>(0.0));
  const std::string checkpoint_out_s = cfg_.get_string("general", "checkpoint_out", std::optional<std::string>(""));
  const std::string resume_from_s = cfg_.get_string("general", "resume_from", std::optional<std::string>(""));

  // Run-level scale (P1): optional.
  ModelScale model_scale = ModelScale::Auto;
  if (cfg_.has_section("model") && cfg_.has_key("model", "scale")) {
    model_scale = parse_model_scale(cfg_.get_string("model", "scale"));
  }

  const fs::path input_path = resolve_path(cfg_dir, input_s);
  const fs::path output_dir_general = resolve_path(cfg_dir, outdir_s);

  std::signal(SIGINT, handle_signal);
  std::signal(SIGTERM, handle_signal);
  if (!validate_only) {
    fs::create_directories(output_dir_general);
  } else {
    // Validation should avoid side-effects. We still resolve paths and will
    // report where outputs would be written.
  }

  // results.json path (relative to general output_dir unless absolute)
  fs::path results_json_path = fs::path(results_json_name);
  if (!results_json_path.is_absolute()) {
    results_json_path = (output_dir_general / results_json_path).lexically_normal();
  }

  // checkpoint paths (relative to general output_dir unless absolute)
  fs::path checkpoint_path;
  if (!checkpoint_out_s.empty()) {
    checkpoint_path = fs::path(checkpoint_out_s);
    if (!checkpoint_path.is_absolute()) checkpoint_path = (output_dir_general / checkpoint_path).lexically_normal();
  }

  fs::path resume_from_path;
  if (!resume_from_s.empty()) {
    resume_from_path = fs::path(resume_from_s);
    if (!resume_from_path.is_absolute()) resume_from_path = (output_dir_general / resume_from_path).lexically_normal();
  }

  const std::uint64_t config_hash = fnv1a64_file(cfg_path.string());

  // --- Reader ---
  LammpsDumpReader reader(input_path.string());
  reader.set_follow(follow, follow_poll_ms);
  reader.set_allow_header_change(allow_header_change);
  reader.set_stop_flag(&g_stop);

  // Profiling accumulators
  double t_reader = 0.0;
  double t_group_setup = 0.0;

  // Read first frame (needed for canonical id mapping + static group setup)
  Frame first_frame;
  {
    ScopedTimer tt(&t_reader);
    const bool ok = reader.next(first_frame);
    if (!ok) throw std::runtime_error("input contains 0 complete frames");
  }

  // --- Groups/Selections ---
  std::unique_ptr<GroupRegistry> group_registry;
  {
    ScopedTimer tg(&t_group_setup);
    std::unordered_map<std::string, std::string> defs;
    if (cfg_.has_section("groups")) {
      const auto& sec = cfg_.section("groups");
      defs.reserve(sec.size());
      for (const auto& [k, v] : sec) defs.emplace(k, v);
    }
    group_registry = std::make_unique<GroupRegistry>(std::move(defs));
    group_registry->prepare_static(first_frame);
  }

  // --- Mapping spec (K6): parsed early for capability validation and checkpoint meta.
  pilots::alg::mapping::MappingSpec mapping_spec;
  mapping_spec.model_scale = model_scale_name(model_scale);
  mapping_spec.mode = pilots::alg::mapping::parse_mapping_mode(
      cfg_.get_string("mapping", "mode", std::optional<std::string>("none")));
  mapping_spec.position = pilots::alg::mapping::parse_mapping_position(
      cfg_.get_string("mapping", "position", std::optional<std::string>("com_geom")));
  mapping_spec.source_group = cfg_.get_string("mapping", "source_group", std::optional<std::string>(""));

  if (!mapping_spec.source_group.empty()) {
    if (!group_registry->has(mapping_spec.source_group)) {
      throw std::runtime_error("mapping.source_group refers to unknown group: '" + mapping_spec.source_group + "'");
    }
    if (group_registry->is_dynamic(mapping_spec.source_group)) {
      throw std::runtime_error(
          "mapping.source_group must be static (identity-consistent), but group '" + mapping_spec.source_group + "' is dynamic");
    }
  }

  // mapping.file path is only required when mode=file.
  if (mapping_spec.mode == pilots::alg::mapping::MappingMode::File) {
    const std::string mp = cfg_.get_string("mapping.file", "path", std::optional<std::string>(""));
    if (mp.empty()) {
      throw std::runtime_error("[mapping] mode=file requires [mapping.file] path=...");
    }
    mapping_spec.file_path = resolve_path(cfg_dir, mp);
    if (!fs::exists(mapping_spec.file_path)) {
      throw std::runtime_error("mapping.file.path does not exist: " + mapping_spec.file_path.string());
    }
    // Include mapping file content hash in the spec hash for reproducibility.
    mapping_spec.file_hash_fnv1a64_hex = hex_u64(pilots::fnv1a64_file(mapping_spec.file_path.string()));
  }

  const bool mapping_enabled = mapping_spec.enabled();
  const std::string mapping_spec_hash_hex = mapping_enabled ? hex_u64(mapping_spec.spec_hash_fnv1a64()) : std::string();

  // Built later (after topology & mol indices).
  std::optional<pilots::alg::mapping::BeadMapping> bead_mapping;
  std::optional<pilots::alg::graph::EdgeList> bead_graph_edges;
  std::optional<pilots::output::PolymerClassifierAudit> polymer_classifier_audit;

  // --- Scan + plan measures (P1-B2/B3) ---
  struct MeasureInstance {
    std::string section;
    std::string instance;
    std::string type;
    MeasureCapabilities caps;
  };

  MeasureBuildEnv env;
  env.cfg_dir = cfg_dir;
  env.input_path = input_path;
  env.output_dir_general = output_dir_general;
  env.dt = dt;
  env.follow = follow;
  env.group_registry = group_registry.get();
  env.first_frame = &first_frame;
  env.model_scale = model_scale;
  env.mapping_available = mapping_enabled;
  env.dry_run = validate_only;

  std::vector<MeasureInstance> planned;
  {
    const auto secs = cfg_.section_names();
    for (const auto& sec : secs) {
      if (!starts_with(sec, "measure.")) continue;
      const std::string instance = sec.substr(std::string("measure.").size());
      if (instance.empty()) {
        throw std::runtime_error("invalid measure section name: [" + sec + "]");
      }

      const bool enabled = cfg_.get_bool(sec, "enabled", std::optional<bool>(true));
      if (!enabled) continue;

      const std::string type = cfg_.get_string(sec, "type", std::optional<std::string>(instance));

      const auto& factory = MeasureRegistry::instance().require(type);
      MeasureCapabilities caps = factory.caps(cfg_, sec, instance, env);

      // --- Runner-side validation (B3) ---
      if (!caps.scale.supports(model_scale)) {
        throw std::runtime_error("measure '" + instance + "' (type='" + type + "') does not support model.scale='" + model_scale_name(model_scale) + "'");
      }

      if (caps.mapping.requires_mapping && !env.mapping_available) {
        throw std::runtime_error(
            "measure '" + instance + "' (type='" + type + "') requires mapping, but mapping is disabled. "
            "Enable it by setting [mapping] mode=identity|by_mol|file");
      }

      // Selection policy validation is centralized later, after SelectionProvider
      // is constructed. This avoids partial checks that only work for named groups
      // (and would miss dynamic expressions like "A & dynamicB").

      // Field dependency checks on the first frame.
      for (const auto& f : caps.requires_dfields) {
        try {
          (void)first_frame.require_dfield(f);
        } catch (const std::exception& e) {
          throw std::runtime_error("measure '" + instance + "' (type='" + type + "') missing required dfield '" + f + "': " + e.what());
        }
      }
      for (const auto& f : caps.requires_i64fields) {
        try {
          (void)first_frame.require_i64field(f);
        } catch (const std::exception& e) {
          throw std::runtime_error("measure '" + instance + "' (type='" + type + "') missing required i64field '" + f + "': " + e.what());
        }
      }
      for (const auto& f : caps.requires_intfields) {
        try {
          (void)first_frame.require_intfield(f);
        } catch (const std::exception& e) {
          throw std::runtime_error("measure '" + instance + "' (type='" + type + "') missing required intfield '" + f + "': " + e.what());
        }
      }

      planned.push_back(MeasureInstance{sec, instance, type, std::move(caps)});
    }
  }

  if (planned.empty()) {
    std::cerr << "[PILOTS] no enabled measures; nothing to do.\n";
    return 0;
  }

  // Selection dependencies (C): if any measure uses topo_group != all, we must have bonds loaded.
  bool selection_needs_topology_bonds = false;
  for (const auto& mi : planned) {
    const std::string topo_group_ref = cfg_.get_string(mi.section, "topo_group", std::optional<std::string>("all"));
    std::string t = topo_group_ref;
    for (auto& c : t) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
    if (!t.empty() && t != "all") {
      if (!cfg_.has_section("topo_groups")) {
        throw std::runtime_error(
            "measure '" + mi.instance + "' (type='" + mi.type + "') requests topo_group='" + topo_group_ref +
            "' but no [topo_groups] section is configured");
      }
      selection_needs_topology_bonds = true;
      break;
    }
  }

  // --- Dependency aggregation: topology sections (B4) ---
  std::vector<std::string> required_topo_sections;
  {
    for (const auto& mi : planned) {
      for (const auto& s : mi.caps.requires_topology_sections) {
        required_topo_sections.push_back(s);
      }
    }

    // Mapping-driven dependencies (K6). Mapping is defined at Runner level, not by measures,
    // but it may require additional topology sections.
    if (mapping_enabled) {
      if (mapping_spec.position == pilots::alg::mapping::MappingPosition::ComMass) {
        if (!first_frame.has_type) {
          throw std::runtime_error("mapping.position=com_mass requires per-atom 'type' to be present in dump (or future topology Atoms parsing)");
        }
        required_topo_sections.push_back("masses");
      }

      if (mapping_spec.mode == pilots::alg::mapping::MappingMode::ByMol && !first_frame.has_mol) {
        const bool derive_mol = cfg_.get_bool("topology", "derive_mol_from_bonds", std::optional<bool>(false));
        if (!derive_mol) {
          throw std::runtime_error(
              "mapping.mode=by_mol requires 'mol' in dump; or set topology.derive_mol_from_bonds=true and provide bonds topology");
        }
        required_topo_sections.push_back("bonds");
      }
    }

    if (selection_needs_topology_bonds) {
      required_topo_sections.push_back("bonds");
    }
    // normalize
    for (auto& s : required_topo_sections) {
      for (auto& c : s) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
    }
    std::sort(required_topo_sections.begin(), required_topo_sections.end());
    required_topo_sections.erase(std::unique(required_topo_sections.begin(), required_topo_sections.end()), required_topo_sections.end());
  }

  // User override: topology.load = ...
  std::vector<std::string> topo_load_sections;
  if (cfg_.has_section("topology") && cfg_.has_key("topology", "load")) {
    topo_load_sections = parse_csv_list(cfg_.get_string("topology", "load"));
  } else {
    topo_load_sections = required_topo_sections;
  }

  TopologyLoadMask topo_mask = to_topology_mask(topo_load_sections);

  if (selection_needs_topology_bonds && !topo_mask.bonds) {
    throw std::runtime_error(
        "at least one measure uses topo_group != all, which requires topology bonds, but topology.load does not include bonds");
  }

  if (mapping_enabled) {
    if (mapping_spec.position == pilots::alg::mapping::MappingPosition::ComMass && !topo_mask.masses) {
      throw std::runtime_error(
          "mapping.position=com_mass requires topology masses, but topology.load does not include masses");
    }
    if (mapping_spec.mode == pilots::alg::mapping::MappingMode::ByMol && !first_frame.has_mol) {
      const bool derive_mol = cfg_.get_bool("topology", "derive_mol_from_bonds", std::optional<bool>(false));
      if (derive_mol && !topo_mask.bonds) {
        throw std::runtime_error(
            "mapping.mode=by_mol with topology.derive_mol_from_bonds=true requires topology bonds, but topology.load does not include bonds");
      }
    }
  }

  // --- Topology loading (optional; only if needed/forced) ---
  // Metadata for results.json / checkpoint meta.
  std::string topology_file_for_index;
  std::string topology_format_for_index;
  std::string topology_hash_for_index;
  std::vector<std::string> topology_loaded_sections_for_index;
  const bool topology_derive_mol_from_bonds = cfg_.has_section("topology")
      ? cfg_.get_bool("topology", "derive_mol_from_bonds", std::optional<bool>(false))
      : false;
  std::optional<Topology> topo;
  const bool need_topology = topo_mask.bonds || topo_mask.masses || topo_mask.angles || topo_mask.dihedrals || topo_mask.impropers;
  if (need_topology) {
    if (!cfg_.has_section("topology") || !cfg_.has_key("topology", "file")) {
      throw std::runtime_error("topology sections requested but no [topology] file is configured");
    }
    const std::string topo_file_s = cfg_.get_string("topology", "file", std::optional<std::string>(""));
    if (topo_file_s.empty()) {
      throw std::runtime_error("topology sections requested but topology.file is empty");
    }

    const std::string topo_format_s = cfg_.get_string("topology", "format", std::optional<std::string>("bond_table"));
    const fs::path topo_path = resolve_path(cfg_dir, topo_file_s);

    std::string fmt = topo_format_s;
    for (auto& c : fmt) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));

    if (fmt == "bond_table" || fmt == "table" || fmt == "bonds") {
      // bond_table can only provide bonds.
      if (topo_mask.masses || topo_mask.angles || topo_mask.dihedrals || topo_mask.impropers) {
        throw std::runtime_error(
            "topology.format=bond_table only supports bonds, but topology.load requested additional sections (masses/angles/dihedrals/impropers)");
      }
      topo = read_bond_table(topo_path);
    } else if (fmt == "lammps_data" || fmt == "lammps" || fmt == "data") {
      topo = read_lammps_data_topology(topo_path, topo_mask);
    } else {
      throw std::runtime_error("unknown topology.format: '" + topo_format_s + "' (use bond_table or lammps_data)");
    }

    topo->finalize_from_frame(first_frame);

    topology_file_for_index = topo_path.string();
    topology_format_for_index = topo_format_s;
    {
      const std::uint64_t h = fnv1a64_file(topo_path.string());
      std::ostringstream ss;
      ss << std::hex << std::setw(16) << std::setfill('0') << h;
      topology_hash_for_index = ss.str();
    }

    topology_loaded_sections_for_index = topo->loaded_sections();
  }

  // --- Molecule index (mol-first, B4) ---
  std::optional<MoleculeIndex> mol_index;
  std::optional<std::vector<std::int64_t>> derived_mol_ids; // for mapping / audit
  if (first_frame.has_mol) {
    mol_index.emplace();
    mol_index->build_from_frame(first_frame);
  } else {
    if (topology_derive_mol_from_bonds) {
      if (!topo || !topo->has_bonds()) {
        throw std::runtime_error("topology.derive_mol_from_bonds=true requires topology bonds, but no bonds were loaded");
      }
      derived_mol_ids = topo->derived_mol_ids_from_bonds(first_frame);
      mol_index.emplace();
      mol_index->build_from_mol_vector(*derived_mol_ids, first_frame.natoms);
    }
  }

  SystemContext sysctx;
  sysctx.topology = topo ? &(*topo) : nullptr;
  sysctx.molecules = mol_index ? &(*mol_index) : nullptr;

  // --- BondGraph (D3): by default, mirror all_bonds adjacency when bonds are loaded.
  std::optional<BondGraph> bond_graph;
  if (sysctx.has_topology() && sysctx.topology->has_bonds()) {
    std::string bond_graph_mode = cfg_.get_string("topology", "bond_graph", std::optional<std::string>("all_bonds"));
    for (auto& c : bond_graph_mode) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
    if (!(bond_graph_mode == "all_bonds" || bond_graph_mode == "all")) {
      throw std::runtime_error("unsupported topology.bond_graph='" + bond_graph_mode + "' (supported: all_bonds)");
    }

    std::vector<int> keep_bond_types;
    if (cfg_.has_key("topology", "bond_graph.bond_types")) {
      const auto s = cfg_.get_list("topology", "bond_graph.bond_types", std::optional<std::string>(""));
      keep_bond_types.reserve(s.size());
      for (const auto& x : s) {
        if (x.empty()) continue;
        keep_bond_types.push_back(std::stoi(x));
      }
    }

    std::string atom_group_filter = cfg_.get_string("topology", "bond_graph.atom_group", std::optional<std::string>(""));
    if (atom_group_filter == "all") atom_group_filter.clear();

    std::vector<unsigned char> keep_atom;
    if (!atom_group_filter.empty()) {
      if (!env.group_registry->has(atom_group_filter)) {
        throw std::runtime_error("topology.bond_graph.atom_group refers to unknown atom group '" + atom_group_filter + "'");
      }
      if (env.group_registry->is_dynamic(atom_group_filter)) {
        throw std::runtime_error("topology.bond_graph.atom_group must be static, but group '" + atom_group_filter + "' is dynamic");
      }
      keep_atom.assign(first_frame.natoms, 0);
      const auto v = env.group_registry->get_view(atom_group_filter, first_frame, 0);
      for (std::size_t i : v.idx) {
        if (i < keep_atom.size()) keep_atom[i] = 1;
      }
    }

    bond_graph.emplace();
    bond_graph->mode = bond_graph_mode;
    bond_graph->natoms = first_frame.natoms;
    bond_graph->bond_types_filter = keep_bond_types;
    bond_graph->atom_group_filter = atom_group_filter;
    bond_graph->adjacency.assign(first_frame.natoms, {});

    auto type_ok = [&](int t) -> bool {
      if (keep_bond_types.empty()) return true;
      for (int kk : keep_bond_types) {
        if (kk == t) return true;
      }
      return false;
    };

    for (const auto& b : sysctx.topology->bonds) {
      if (!type_ok(b.type)) continue;
      if (!keep_atom.empty()) {
        if (b.i >= keep_atom.size() || b.j >= keep_atom.size()) continue;
        if (!(keep_atom[b.i] && keep_atom[b.j])) continue;
      }
      bond_graph->adjacency[b.i].push_back(b.j);
      bond_graph->adjacency[b.j].push_back(b.i);
    }
    // Optional cleanup: sort/unique adjacency lists for determinism.
    for (auto& nbrs : bond_graph->adjacency) {
      std::sort(nbrs.begin(), nbrs.end());
      nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
    }

    sysctx.bond_graph = &(*bond_graph);
  }

  // --- Bead mapping (K6) ---
  if (mapping_enabled) {
    std::vector<std::size_t> atoms_filter;
    if (!mapping_spec.source_group.empty()) {
      const auto v = group_registry->get_view(mapping_spec.source_group, first_frame, /*frame_idx*/0);
      atoms_filter.assign(v.idx.begin(), v.idx.end());
    }

    const std::span<const std::size_t> filter_span(atoms_filter.data(), atoms_filter.size());

    switch (mapping_spec.mode) {
      case pilots::alg::mapping::MappingMode::Identity:
        bead_mapping = pilots::alg::mapping::BeadMapping::build_identity(first_frame, mapping_spec, filter_span);
        break;
      case pilots::alg::mapping::MappingMode::ByMol:
        if (first_frame.has_mol) {
          bead_mapping = pilots::alg::mapping::BeadMapping::build_by_mol(first_frame, mapping_spec, filter_span);
        } else if (derived_mol_ids) {
          const std::span<const std::int64_t> mol_span(derived_mol_ids->data(), derived_mol_ids->size());
          bead_mapping = pilots::alg::mapping::BeadMapping::build_by_mol_ids(first_frame, mol_span, mapping_spec, filter_span);
        } else {
          throw std::runtime_error(
              "mapping.mode=by_mol requires mol ids (dump 'mol' field or topology.derive_mol_from_bonds=true with bonds loaded)");
        }
        break;
      case pilots::alg::mapping::MappingMode::File:
        bead_mapping = pilots::alg::mapping::BeadMapping::build_from_file(first_frame, mapping_spec, mapping_spec.file_path, filter_span);
        break;
      case pilots::alg::mapping::MappingMode::None:
        break;
    }

    if (!bead_mapping) {
      throw std::runtime_error("mapping.enabled but mapping definition was not created (internal error)");
    }

    // Derive per-atom/bead masses for audit and com_mass (if requested).
    bead_mapping->compute_masses_from(first_frame, topo ? &(*topo) : nullptr);
    if (mapping_spec.position == pilots::alg::mapping::MappingPosition::ComMass && !bead_mapping->has_masses_for_com()) {
      throw std::runtime_error(
          "mapping.position=com_mass requires topology masses (Masses section loaded) and dump 'type' field");
    }

    sysctx.mapping = &(*bead_mapping);
  }

  // --- Polymer classifier (K7) ---
  {
    const bool pc_enabled = cfg_.get_bool("polymer_classifier", "enabled", std::optional<bool>(true));
    if (pc_enabled) {
      auto lower_ascii = [](std::string s) {
        for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return s;
      };

      std::string src = cfg_.get_string("polymer_classifier", "graph_source", std::optional<std::string>("auto"));
      src = lower_ascii(src);
      if (src.empty()) src = "auto";

      // Determine the actual graph source.
      std::string graph_source = src;
      if (graph_source == "auto") {
        if (bead_mapping && topo && topo->has_bonds()) {
          graph_source = "bead_graph";
        } else if (bond_graph) {
          graph_source = "bond_graph";
        } else if (topo && topo->has_bonds()) {
          graph_source = "topology_graph";
        } else {
          graph_source = "none";
        }
      }

      pilots::output::PolymerClassifierAudit audit;
      audit.graph_source = graph_source;

      if (graph_source == "none") {
        audit.tags = {"unavailable"};
        audit.hint = "polymer classifier skipped: no topology bonds were loaded (topology.load must include bonds)";
        audit.n_components = 0;
        audit.largest_component_size = 0;
      } else {
        std::optional<pilots::alg::graph::GraphView> gv;
        if (graph_source == "topology_graph") {
          if (!topo || !topo->has_bonds()) {
            throw std::runtime_error("polymer_classifier.graph_source=topology_graph requires topology bonds, but no bonds were loaded");
          }
          gv.emplace(topo->adjacency);
        } else if (graph_source == "bond_graph") {
          if (!bond_graph) {
            throw std::runtime_error("polymer_classifier.graph_source=bond_graph requires bonds, but bond_graph is not available");
          }
          gv.emplace(bond_graph->adjacency);
        } else if (graph_source == "bead_graph") {
          if (!bead_mapping) {
            throw std::runtime_error("polymer_classifier.graph_source=bead_graph requires mapping, but mapping is disabled");
          }
          if (!topo || !topo->has_bonds()) {
            throw std::runtime_error("polymer_classifier.graph_source=bead_graph requires topology bonds, but no bonds were loaded");
          }
          bead_graph_edges = pilots::alg::mapping::build_bead_graph_from_bonds(*topo, *bead_mapping);
          gv.emplace(*bead_graph_edges);
        } else {
          throw std::runtime_error("invalid polymer_classifier.graph_source='" + src + "' (expected auto|topology_graph|bond_graph|bead_graph)");
        }

        const auto res = pilots::alg::polymer::classify_graph(*gv, graph_source, model_scale_name(model_scale));
        audit.tags = res.tags;
        audit.hint = res.hint;
        audit.n_components = res.n_components();
        audit.largest_component_size = res.largest_component_size();
        audit.components.reserve(res.components.size());
        for (const auto& c : res.components) {
          pilots::output::PolymerClassifierComponentAudit cc;
          cc.id = c.id;
          cc.size = c.size;
          cc.n_ends = c.n_ends;
          cc.n_branch_points = c.n_branch_points;
          cc.cycle_rank = c.cycle_rank;
          cc.tag = c.tag;
          audit.components.push_back(std::move(cc));
        }
      }

      polymer_classifier_audit = std::move(audit);
    }
  }

  // --- TopoGroups (C) ---
  std::optional<TopoGroupRegistry> topo_groups;
  if (cfg_.has_section("topo_groups")) {
    topo_groups.emplace(cfg_.section("topo_groups"));
    if (sysctx.has_topology() && sysctx.topology->has_bonds()) {
      topo_groups->prepare_static(first_frame, *sysctx.topology, *env.group_registry);
    }
  }

  // Unified selection provider (AtomGroup + TopoGroup + A/T combine).
  SelectionProvider selection_provider(env.group_registry, topo_groups ? &(*topo_groups) : nullptr,
                                       sysctx.topology);
  env.selection_provider = &selection_provider;

  // --- Centralized selection-policy validation (B3/C2) ---
  // Measures can still perform internal checks, but the Runner should fail-fast before
  // instantiation when a measure declares it requires static selections.
  for (const auto& mi : planned) {
    if (mi.caps.selection_policy != SelectionPolicy::RequireStatic) continue;

    // Standard per-measure selection keys (group/topo_group). These are part of the platform
    // contract: any measure may use them for its "main" selection.
    {
      const std::string group_ref = cfg_.get_string(mi.section, "group", std::optional<std::string>("all"));
      const std::string topo_group_ref = cfg_.get_string(mi.section, "topo_group", std::optional<std::string>("all"));
      if (selection_provider.is_dynamic_spec(group_ref, topo_group_ref)) {
        throw std::runtime_error(
            "measure '" + mi.instance + "' (type='" + mi.type + "') requires a static selection, but the combined selection depends on a dynamic group/topo_group");
      }
    }

    // Additional group references declared by the factory (e.g., drift_group for MSD).
    for (const auto& gref : mi.caps.group_refs) {
      if (selection_provider.is_dynamic_spec(gref, "all")) {
        throw std::runtime_error(
            "measure '" + mi.instance + "' (type='" + mi.type + "') requires a static selection, but group '" + gref + "' is dynamic (or uses a dynamic selector)");
      }
    }
  }

  // Validate topology availability per measure now that we know what is loaded.
  for (const auto& mi : planned) {
    for (auto s : mi.caps.requires_topology_sections) {
      for (auto& c : s) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
      if (s == "bonds" || s == "bond") {
        if (!sysctx.has_topology() || !sysctx.topology->has_bonds()) {
          throw std::runtime_error("measure '" + mi.instance + "' (type='" + mi.type + "') requires topology section 'bonds', but bonds were not loaded");
        }
      } else if (s == "masses" || s == "mass") {
        if (!sysctx.has_topology() || !sysctx.topology->has_masses()) {
          throw std::runtime_error("measure '" + mi.instance + "' (type='" + mi.type + "') requires topology section 'masses', but masses were not loaded");
        }
      } else if (s == "angles" || s == "angle") {
        if (!sysctx.has_topology() || !sysctx.topology->has_angles()) {
          throw std::runtime_error("measure '" + mi.instance + "' (type='" + mi.type + "') requires topology section 'angles', but angles were not loaded");
        }
      } else if (s == "dihedrals" || s == "dihedral") {
        if (!sysctx.has_topology() || !sysctx.topology->has_dihedrals()) {
          throw std::runtime_error("measure '" + mi.instance + "' (type='" + mi.type + "') requires topology section 'dihedrals', but dihedrals were not loaded");
        }
      } else if (s == "impropers" || s == "improper") {
        if (!sysctx.has_topology() || !sysctx.topology->has_impropers()) {
          throw std::runtime_error("measure '" + mi.instance + "' (type='" + mi.type + "') requires topology section 'impropers', but impropers were not loaded");
        }
      } else {
        throw std::runtime_error("measure '" + mi.instance + "' (type='" + mi.type + "') requires unknown topology section '" + s + "'");
      }
    }
  }

  // --- Instantiate measures (B2) ---
  std::vector<std::unique_ptr<IMeasure>> measures;
  std::vector<output::MeasureProfiling> mp;
  measures.reserve(planned.size());
  mp.assign(planned.size(), output::MeasureProfiling{});

  for (const auto& mi : planned) {
    const auto& factory = MeasureRegistry::instance().require(mi.type);
    measures.emplace_back(factory.create(cfg_, mi.section, mi.instance, env, sysctx));
  }

  for (std::size_t i = 0; i < measures.size(); ++i) {
    ScopedTimer tm(&mp[i].on_start_s);
    measures[i]->on_start(first_frame, sysctx);
  }

  if (validate_only) {
    std::cerr << "[PILOTS] validation OK (no trajectory processing performed)\n"
              << "         input=" << input_path.string() << "\n"
              << "         measures=" << measures.size() << " output_dir=" << output_dir_general.string() << "\n";
    if (sysctx.has_topology()) {
      std::cerr << "         topology_sections=";
      const auto secs = sysctx.topology->loaded_sections();
      for (std::size_t i = 0; i < secs.size(); ++i) {
        if (i) std::cerr << ",";
        std::cerr << secs[i];
      }
      std::cerr << "\n";
    }
    if (mapping_enabled) {
      std::cerr << "         mapping=enabled mode=" << pilots::alg::mapping::mapping_mode_name(mapping_spec.mode)
                << " position=" << pilots::alg::mapping::mapping_position_name(mapping_spec.position) << "\n";
    }
    return 0;
  }

  // --- Resume (optional) ---
  bool resumed = false;
  double wall_base = 0.0;
  std::size_t frame_idx = 0;

  if (!resume_from_path.empty() && fs::exists(resume_from_path)) {
    const ResumeMeta meta = load_checkpoint_into(resume_from_path, config_hash, input_path, mapping_spec_hash_hex, measures, mp);

    // If v2+ checkpoint carries topology meta, sanity-check consistency.
    if (!meta.topology_hash_fnv1a64_hex.empty()) {
      if (meta.topology_hash_fnv1a64_hex != topology_hash_for_index) {
        throw std::runtime_error("checkpoint topology hash mismatch: checkpoint=" + meta.topology_hash_fnv1a64_hex + " current=" + topology_hash_for_index);
      }
      {
        auto a = meta.topology_loaded_sections;
        auto b = topology_loaded_sections_for_index;
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        if (a != b) {
          throw std::runtime_error("checkpoint topology.loaded_sections mismatch");
        }
      }
      if (meta.topology_derive_mol_from_bonds != topology_derive_mol_from_bonds) {
        throw std::runtime_error("checkpoint topology.derive_mol_from_bonds mismatch");
      }
    }

    wall_base = meta.wall_seconds;
    t_reader += meta.reader_seconds;
    t_group_setup += meta.group_setup_seconds;

    reader.seek_offset(meta.reader_offset);
    frame_idx = static_cast<std::size_t>(meta.frame_idx);
    resumed = true;

    std::cerr << "[PILOTS] resumed from checkpoint: " << resume_from_path.string() << "\n"
              << "         next_frame_index=" << frame_idx << " reader_offset=" << meta.reader_offset << "\n";
  }

  // --- Results index scaffolding ---
  output::ResultsIndex idx;
  idx.pilots_version = PILOTS_VERSION_STR;
  idx.schema_version = output::RESULTS_SCHEMA_VERSION;
  idx.model_scale = model_scale_name(model_scale);
  idx.dt_unit = "simulation_time";
  idx.config_path = cfg_path.string();
  idx.config_hash_fnv1a64_hex = hex_u64(config_hash);
  idx.input_path = input_path.string();
  // File fingerprints (F): for reproducibility, record size/mtime and (when stable) a content hash.
  // - config: always hashed (small, immutable)
  // - topology: hashed when present
  // - input dump: hashed by default when NOT follow; in follow mode the file may still be growing
  //              and hashing can be misleading/expensive.
  const bool hash_input = cfg_.get_bool("general", "hash_input", std::optional<bool>(!follow));
  idx.config_fingerprint = make_fingerprint(cfg_path, /*compute_hash=*/true, config_hash);
  idx.input_fingerprint = make_fingerprint(input_path, /*compute_hash=*/hash_input);

  // Mapping audit (K6 + F): records bead definition requests in a stable hash.
  {
    output::MappingAudit ma;
    ma.enabled = mapping_enabled;
    ma.mode = pilots::alg::mapping::mapping_mode_name(mapping_spec.mode);
    ma.position = pilots::alg::mapping::mapping_position_name(mapping_spec.position);
    ma.source_group = mapping_spec.source_group;
    ma.file = mapping_spec.file_path.empty() ? std::string("") : mapping_spec.file_path.string();
    ma.file_hash_fnv1a64_hex = mapping_spec.file_hash_fnv1a64_hex;
    // Only provide a spec hash when mapping is enabled. This avoids confusion
    // in downstream tooling where enabled=false but a non-empty hash looks like
    // an active mapping definition.
    ma.spec_hash_fnv1a64_hex = mapping_enabled ? hex_u64(mapping_spec.spec_hash_fnv1a64()) : std::string();
    idx.mapping = std::move(ma);
  }

  if (polymer_classifier_audit) {
    idx.polymer_classifier = *polymer_classifier_audit;
  }
  idx.output_dir = output_dir_general.string();
  idx.dt = dt;
  idx.follow = follow;
  idx.flush_every_frames = flush_every_frames;
  idx.flush_every_seconds = flush_every_seconds;
  idx.checkpoint_out = checkpoint_path.empty() ? std::string("") : checkpoint_path.string();
  idx.resume_from = resume_from_path.empty() ? std::string("") : resume_from_path.string();

  // Topology metadata (D1/D2)
  idx.topology_file = topology_file_for_index;
  idx.topology_format = topology_format_for_index;
  idx.topology_hash_fnv1a64_hex = topology_hash_for_index;
  if (!topology_file_for_index.empty()) {
    // Avoid hashing twice: we already computed topology_hash_for_index during load.
    idx.topology_fingerprint = make_fingerprint(fs::path(topology_file_for_index), /*compute_hash=*/false);
    if (!topology_hash_for_index.empty()) {
      idx.topology_fingerprint.hash_kind = "fnv1a64";
      idx.topology_fingerprint.hash_fnv1a64_hex = topology_hash_for_index;
      idx.topology_fingerprint.hash_computed = true;
    }
  } else {
    idx.topology_fingerprint = output::FileFingerprint{};
    idx.topology_fingerprint.path = "";
    idx.topology_fingerprint.hash_kind = "none";
    idx.topology_fingerprint.hash_computed = false;
  }
  idx.topology_loaded_sections = topology_loaded_sections_for_index;
  idx.topology_derive_mol_from_bonds = topology_derive_mol_from_bonds;

  // BondGraph audit (D3)
  if (bond_graph) {
    output::BondGraphAudit bga;
    bga.mode = bond_graph->mode;
    bga.atom_group = bond_graph->atom_group_filter;
    bga.bond_types = bond_graph->bond_types_filter;
    std::size_t deg_sum = 0;
    for (const auto& nbrs : bond_graph->adjacency) deg_sum += nbrs.size();
    bga.edges = deg_sum / 2;
    idx.bond_graph = std::move(bga);
  } else {
    idx.bond_graph.reset();
  }

#if PILOTS_HAS_OPENMP
  idx.threads = omp_get_max_threads();
#else
  idx.threads = 1;
#endif

  idx.measures.reserve(measures.size());
  idx.measure_profiling.reserve(measures.size());
  for (std::size_t i = 0; i < measures.size(); ++i) {
    idx.measures.push_back(measures[i]->describe());
    idx.measure_profiling.emplace_back(measures[i]->instance_name(), mp[i]);
  }

  auto to_output_group_audits = [](const auto& infos) {
    std::vector<output::GroupAudit> out;
    out.reserve(infos.size());
    for (const auto& info : infos) {
      output::GroupAudit ga;
      ga.name = info.name;
      ga.expr = info.expr;
      ga.is_dynamic = info.is_dynamic;
      if (!info.is_dynamic) {
        ga.samples = 1;
        ga.size_min = info.static_size;
        ga.size_max = info.static_size;
        ga.size_mean = static_cast<double>(info.static_size);
        ga.size_var = 0.0;
      } else {
        ga.samples = info.size_stats.samples;
        ga.size_min = info.size_stats.min;
        ga.size_max = info.size_stats.max;
        ga.size_mean = info.size_stats.mean;
        ga.size_var = info.size_stats.variance();
      }
      out.push_back(std::move(ga));
    }
    return out;
  };

  auto write_index = [&]() {
    idx.frames_processed = frame_idx;
    idx.wall_seconds = wall_base + total_timer.elapsed_seconds();
    idx.reader_seconds = t_reader;
    idx.group_setup_seconds = t_group_setup;
    idx.reader_offset = reader.tell_offset();

    idx.measure_profiling.clear();
    idx.measure_profiling.reserve(measures.size());
    for (std::size_t i = 0; i < measures.size(); ++i) {
      idx.measure_profiling.emplace_back(measures[i]->instance_name(), mp[i]);
    }
    idx.atom_groups = to_output_group_audits(group_registry->audit());
    idx.topo_groups = topo_groups ? to_output_group_audits(topo_groups->audit()) : std::vector<output::GroupAudit>{};

    // Combined selection audits (F/C3): provide auditable size stats for any combined
    // selection spec that was evaluated (dynamic selections accumulate stats over frames).
    idx.combined_selections.clear();
    {
      const auto sel_audits = selection_provider.audit();
      idx.combined_selections.reserve(sel_audits.size());
      for (const auto& a : sel_audits) {
        output::CombinedSelectionAudit ca;
        ca.key = a.key;
        ca.is_dynamic = a.is_dynamic;
        if (!a.is_dynamic) {
          ca.samples = 1;
          ca.size_min = a.static_size;
          ca.size_max = a.static_size;
          ca.size_mean = static_cast<double>(a.static_size);
          ca.size_var = 0.0;
        } else {
          ca.samples = a.dynamic_size.samples;
          ca.size_min = a.dynamic_size.min;
          ca.size_max = a.dynamic_size.max;
          ca.size_mean = a.dynamic_size.mean;
          ca.size_var = a.dynamic_size.variance();
        }
        idx.combined_selections.push_back(std::move(ca));
      }
    }
    output::write_results_json(results_json_path, idx);
  };

  // Write an initial index snapshot (useful in resume mode)
  write_index();

  // If not resumed, process the already-read first frame as frame 0.
  Frame frame = first_frame;
  if (!resumed) {
    for (std::size_t i = 0; i < measures.size(); ++i) {
      ScopedTimer tm(&mp[i].on_frame_s);
      measures[i]->on_frame(frame, 0);
      mp[i].frames += 1;
    }
    frame_idx = 1;
  }

  // Flush bookkeeping (online/follow)
  std::size_t last_flush_frame = frame_idx;
  double last_flush_t = total_timer.elapsed_seconds();

  auto do_flush = [&]() {
    for (auto& m : measures) {
      m->flush_partial();
    }
    write_index();

    if (!checkpoint_path.empty()) {
      save_checkpoint_from(checkpoint_path,
                           config_hash,
                           input_path,
                           mapping_spec_hash_hex,
                           reader.tell_offset(),
                           static_cast<std::uint64_t>(frame_idx),
                           idx.wall_seconds,
                           t_reader,
                           t_group_setup,
                           topology_file_for_index,
                           topology_format_for_index,
                           topology_hash_for_index,
                           topology_loaded_sections_for_index,
                           topology_derive_mol_from_bonds,
                           measures,
                           mp);
    }
  };

  // --- Main read loop ---
  while (true) {
    if (g_stop) break;

    bool ok = false;
    {
      ScopedTimer tt(&t_reader);
      ok = reader.next(frame);
    }
    if (!ok) break;

    for (std::size_t i = 0; i < measures.size(); ++i) {
      ScopedTimer tm(&mp[i].on_frame_s);
      measures[i]->on_frame(frame, frame_idx);
      mp[i].frames += 1;
    }
    ++frame_idx;

    const bool need_flush_frames = (flush_every_frames > 0) && (frame_idx >= last_flush_frame + flush_every_frames);
    const bool need_flush_time = (flush_every_seconds > 0.0) && ((total_timer.elapsed_seconds() - last_flush_t) >= flush_every_seconds);
    if (need_flush_frames || need_flush_time) {
      do_flush();
      last_flush_frame = frame_idx;
      last_flush_t = total_timer.elapsed_seconds();
    }
  }

  // Finalize measures
  for (std::size_t i = 0; i < measures.size(); ++i) {
    ScopedTimer tm(&mp[i].finalize_s);
    measures[i]->finalize();
  }

  // Final flush/index + final checkpoint
  do_flush();

  if (print_profile) {
    std::cerr << "[PILOTS] profiling\n";
    std::cerr << "  wall_seconds: " << std::setprecision(6) << idx.wall_seconds << "\n";
    std::cerr << "  reader_seconds: " << idx.reader_seconds << "\n";
    std::cerr << "  group_setup_seconds: " << idx.group_setup_seconds << "\n";
    for (std::size_t i = 0; i < measures.size(); ++i) {
      std::cerr << "  measure." << measures[i]->instance_name() << " (type=" << measures[i]->type() << ")"
                << ": on_start=" << mp[i].on_start_s
                << " on_frame=" << mp[i].on_frame_s
                << " finalize=" << mp[i].finalize_s
                << " frames=" << mp[i].frames << "\n";
    }
    std::cerr << "  results_json: " << results_json_path.string() << "\n";
    if (!checkpoint_path.empty()) {
      std::cerr << "  checkpoint: " << checkpoint_path.string() << "\n";
    }
  }

  return 0;
}

} // namespace pilots
