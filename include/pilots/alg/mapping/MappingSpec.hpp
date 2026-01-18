#pragma once

#include <cctype>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>

#include "pilots/util/Hash.hpp"

namespace pilots::alg::mapping {
namespace fs = std::filesystem;

enum class MappingMode {
  None = 0,
  Identity = 1,
  ByMol = 2,
  File = 3,
};

enum class MappingPosition {
  ComMass = 0,
  ComGeom = 1,
  RepresentativeAtom = 2,
};

inline std::string to_lower(std::string s) {
  for (auto& c : s) c = static_cast<char>(::tolower(static_cast<unsigned char>(c)));
  return s;
}

inline MappingMode parse_mapping_mode(std::string s) {
  s = to_lower(std::move(s));
  if (s.empty() || s == "none") return MappingMode::None;
  if (s == "identity") return MappingMode::Identity;
  if (s == "by_mol" || s == "bymol" || s == "mol") return MappingMode::ByMol;
  if (s == "file") return MappingMode::File;
  throw std::runtime_error("invalid mapping.mode: '" + s + "' (use none|identity|by_mol|file)");
}

inline std::string mapping_mode_name(MappingMode m) {
  switch (m) {
    case MappingMode::None: return "none";
    case MappingMode::Identity: return "identity";
    case MappingMode::ByMol: return "by_mol";
    case MappingMode::File: return "file";
  }
  return "none";
}

inline MappingPosition parse_mapping_position(std::string s) {
  s = to_lower(std::move(s));
  if (s.empty() || s == "com_geom" || s == "comgeom") return MappingPosition::ComGeom;
  if (s == "com_mass" || s == "commass") return MappingPosition::ComMass;
  if (s == "representative_atom" || s == "rep_atom" || s == "rep") return MappingPosition::RepresentativeAtom;
  throw std::runtime_error("invalid mapping.position: '" + s + "' (use com_mass|com_geom|representative_atom)");
}

inline std::string mapping_position_name(MappingPosition p) {
  switch (p) {
    case MappingPosition::ComMass: return "com_mass";
    case MappingPosition::ComGeom: return "com_geom";
    case MappingPosition::RepresentativeAtom: return "representative_atom";
  }
  return "com_geom";
}

// Mapping spec is an auditable, serializable definition of how beads are built.
//
// Important: spec_hash_fnv1a64 is computed from canonical_string(), which is a stable
// serialization of all mapping-relevant fields (including mapping file hash when mode=file).
struct MappingSpec {
  // Run-level scale (aa/ua/cg/auto) for audit/reproducibility. Included in spec_hash.
  std::string model_scale = "auto";

  MappingMode mode = MappingMode::None;
  MappingPosition position = MappingPosition::ComGeom;

  // Optional: restrict mapping to a static atom-group (evaluated on the first frame).
  std::string source_group;

  // For mode=file
  fs::path file_path;            // resolved absolute or cfg-relative path
  std::string file_hash_fnv1a64_hex; // empty if not applicable

  bool enabled() const { return mode != MappingMode::None; }

  std::string canonical_string() const {
    // Stable serialization (order fixed). All strings are written as-is; callers should
    // canonicalize paths before setting file_path.
    std::string out;
    out.reserve(256);
    out += "model_scale=" + model_scale;
    out += ";mode=" + mapping_mode_name(mode);
    out += ";position=" + mapping_position_name(position);
    out += ";source_group=" + source_group;
    // IMPORTANT: file_path must NOT be part of the hash identity. The same mapping file
    // content may live under different paths (e.g., different HPC scratch dirs) but should
    // still be considered the same bead definition. The file content hash provides the
    // reproducibility guarantee.
    out += ";file_hash_fnv1a64=" + file_hash_fnv1a64_hex;
    return out;
  }

  std::uint64_t spec_hash_fnv1a64() const {
    const std::string s = canonical_string();
    return pilots::fnv1a64_bytes(s.data(), s.size());
  }
};

} // namespace pilots::alg::mapping
