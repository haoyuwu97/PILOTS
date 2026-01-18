#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pilots/config/IniConfig.hpp"
#include "pilots/core/Frame.hpp"
#include "pilots/core/SystemContext.hpp"
#include "pilots/measures/IMeasure.hpp"
#include "pilots/select/GroupRegistry.hpp"
#include "pilots/select/SelectionProvider.hpp"

namespace pilots {

// --- Model scale (trajectory interpretation) ---
//
// P1 baseline:
// - The platform records a run-level scale declaration (AA/UA/CG/auto).
// - Measures can declare which scales they support and/or require mapping.
// - Mapping is not implemented in this package yet; the capability system is
//   added now to make future extensions clean and scientifically safe.

enum class ModelScale {
  Auto = 0,
  AA = 1,
  UA = 2,
  CG = 3,
};

inline ModelScale parse_model_scale(std::string s) {
  for (auto& c : s) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
  if (s.empty() || s == "auto") return ModelScale::Auto;
  if (s == "aa" || s == "all-atom" || s == "all_atom") return ModelScale::AA;
  if (s == "ua" || s == "united-atom" || s == "united_atom") return ModelScale::UA;
  if (s == "cg" || s == "coarse-grained" || s == "coarse_grained") return ModelScale::CG;
  throw std::runtime_error("invalid model scale: '" + s + "' (use auto|AA|UA|CG)");
}

inline std::string model_scale_name(ModelScale s) {
  switch (s) {
    case ModelScale::Auto: return "auto";
    case ModelScale::AA: return "AA";
    case ModelScale::UA: return "UA";
    case ModelScale::CG: return "CG";
  }
  return "auto";
}

struct ScaleCompatibility {
  bool aa = true;
  bool ua = true;
  bool cg = true;

  bool supports(ModelScale s) const {
    switch (s) {
      case ModelScale::Auto: return true;
      case ModelScale::AA: return aa;
      case ModelScale::UA: return ua;
      case ModelScale::CG: return cg;
    }
    return true;
  }
};

// --- Selection + dependency capabilities ---

enum class SelectionPolicy {
  AllowDynamic = 0,
  RequireStatic = 1,
};

struct MappingRequirements {
  // Placeholder for the upcoming mapping subsystem.
  bool requires_mapping = false;
  bool requires_beads = false;
  bool requires_chain_pos = false;
  bool requires_bead_graph = false;
};

struct MeasureCapabilities {
  SelectionPolicy selection_policy = SelectionPolicy::AllowDynamic;

  // Time-correlation measures typically require consistent particle identity
  // across frames (i.e. canonical ordering by id).
  bool requires_identity_consistent = false;

  // Field dependencies.
  std::vector<std::string> requires_dfields;   // require_dfield(...)
  std::vector<std::string> requires_i64fields; // require_i64field(...)
  std::vector<std::string> requires_intfields; // require_intfield(...)

  // Topology sections needed (names: masses,bonds,angles,dihedrals,impropers,...)
  std::vector<std::string> requires_topology_sections;

  // Run-level scale compatibility.
  ScaleCompatibility scale;

  // Mapping requirements.
  MappingRequirements mapping;

  // Group references (raw strings) used by this measure instance; Runner uses
  // this list to validate static/dynamic policy.
  std::vector<std::string> group_refs;
};

// Context used by measure factories.
struct MeasureBuildEnv {
  std::filesystem::path cfg_dir;
  std::filesystem::path input_path;
  std::filesystem::path output_dir_general;

  double dt = 0.0;
  bool follow = false;

  GroupRegistry* group_registry = nullptr; // already prepared on first frame
  const Frame* first_frame = nullptr;

  // Optional: present when topology/topo_groups are available.
  // Measures that want unified AtomGroup/TopoGroup selection should use this.
  SelectionProvider* selection_provider = nullptr;

  ModelScale model_scale = ModelScale::Auto;
  bool mapping_available = false;

  // When true, factories should avoid side-effects such as creating directories
  // or opening output files. Used for CLI --validate-config.
  bool dry_run = false;
};

struct MeasureFactoryEntry {
  using CapsFn = MeasureCapabilities (*)(const IniConfig&, const std::string& section,
                                         const std::string& instance,
                                         const MeasureBuildEnv& env);

  using CreateFn = std::unique_ptr<IMeasure> (*)(const IniConfig&, const std::string& section,
                                                 const std::string& instance,
                                                 const MeasureBuildEnv& env,
                                                 const SystemContext& sysctx);

  std::string type;
  CapsFn caps = nullptr;
  CreateFn create = nullptr;
};

class MeasureRegistry {
public:
  static MeasureRegistry& instance() {
    static MeasureRegistry r;
    return r;
  }

  void register_factory(MeasureFactoryEntry e) {
    if (e.type.empty()) throw std::runtime_error("MeasureRegistry: factory type is empty");
    if (!e.caps || !e.create) throw std::runtime_error("MeasureRegistry: factory entry missing function pointers for type='" + e.type + "'");
    auto it = factories_.find(e.type);
    if (it != factories_.end()) {
      throw std::runtime_error("MeasureRegistry: duplicate factory registration for type='" + e.type + "'");
    }
    factories_.emplace(e.type, std::move(e));
  }

  bool has(const std::string& type) const {
    return factories_.find(type) != factories_.end();
  }

  const MeasureFactoryEntry& require(const std::string& type) const {
    auto it = factories_.find(type);
    if (it == factories_.end()) {
      std::vector<std::string> keys;
      keys.reserve(factories_.size());
      for (const auto& kv : factories_) keys.push_back(kv.first);
      std::sort(keys.begin(), keys.end());
      std::string msg = "MeasureRegistry: unknown measure type '" + type + "'. Registered types:";
      for (const auto& k : keys) msg += " " + k;
      throw std::runtime_error(msg);
    }
    return it->second;
  }

  std::vector<std::string> registered_types() const {
    std::vector<std::string> out;
    out.reserve(factories_.size());
    for (const auto& kv : factories_) out.push_back(kv.first);
    std::sort(out.begin(), out.end());
    return out;
  }

private:
  std::unordered_map<std::string, MeasureFactoryEntry> factories_;
};

// Helper to register a measure factory from a single translation unit.
//
// Intended usage (inside one .cpp file per measure):
//
//   static MeasureRegistrar reg("msd", &caps_fn, &create_fn);
//
// This enables the "add one .cpp file" workflow when CMake compiles
// all sources under src/measures/.
class MeasureRegistrar {
public:
  MeasureRegistrar(const char* type,
                   MeasureFactoryEntry::CapsFn caps,
                   MeasureFactoryEntry::CreateFn create) {
    MeasureFactoryEntry e;
    e.type = std::string(type ? type : "");
    e.caps = caps;
    e.create = create;
    MeasureRegistry::instance().register_factory(std::move(e));
  }
};

inline bool starts_with(std::string_view s, std::string_view pref) {
  return s.size() >= pref.size() && s.substr(0, pref.size()) == pref;
}

} // namespace pilots
