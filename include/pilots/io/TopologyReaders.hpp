#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "pilots/topology/Topology.hpp"

namespace pilots {
namespace fs = std::filesystem;

// Which sections to load from a topology file.
//
// For LAMMPS data files:
//   masses, bonds, angles, dihedrals, impropers
// For bond tables:
//   bonds only (other sections are unsupported).
struct TopologySectionMask {
  bool masses = false;
  bool bonds = false;
  bool angles = false;
  bool dihedrals = false;
  bool impropers = false;
};

namespace detail_topo {
inline std::string trim(std::string_view in) {
  std::size_t b = 0;
  while (b < in.size() && (in[b] == ' ' || in[b] == '\t' || in[b] == '\r' || in[b] == '\n')) ++b;
  std::size_t e = in.size();
  while (e > b && (in[e - 1] == ' ' || in[e - 1] == '\t' || in[e - 1] == '\r' || in[e - 1] == '\n')) --e;
  return std::string(in.substr(b, e - b));
}

inline std::string strip_comment(std::string_view s) {
  const auto pos = s.find('#');
  if (pos == std::string_view::npos) return std::string(s);
  return std::string(s.substr(0, pos));
}

inline std::vector<std::string> split_ws(std::string_view s) {
  std::vector<std::string> out;
  std::size_t i = 0;
  while (i < s.size()) {
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t' || s[i] == '\r' || s[i] == '\n')) ++i;
    if (i >= s.size()) break;
    std::size_t j = i;
    while (j < s.size() && !(s[j] == ' ' || s[j] == '\t' || s[j] == '\r' || s[j] == '\n')) ++j;
    out.emplace_back(s.substr(i, j - i));
    i = j;
  }
  return out;
}

inline bool starts_with_alpha(std::string_view s) {
  for (char c : s) {
    if (c == ' ' || c == '\t' || c == '\r' || c == '\n') continue;
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
  }
  return false;
}

inline std::int64_t to_i64(const std::string& x, const std::string& what) {
  try {
    std::size_t pos = 0;
    long long v = std::stoll(x, &pos);
    if (pos != x.size()) throw std::invalid_argument("trailing");
    return static_cast<std::int64_t>(v);
  } catch (...) {
    throw std::runtime_error("Topology: invalid integer token for " + what + ": '" + x + "'");
  }
}

inline int to_int(const std::string& x, const std::string& what) {
  const auto v = to_i64(x, what);
  if (v < static_cast<std::int64_t>(std::numeric_limits<int>::min()) ||
      v > static_cast<std::int64_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("Topology: integer overflow for " + what + ": '" + x + "'");
  }
  return static_cast<int>(v);
}

inline double to_f64(const std::string& x, const std::string& what) {
  try {
    std::size_t pos = 0;
    double v = std::stod(x, &pos);
    if (pos != x.size()) throw std::invalid_argument("trailing");
    return v;
  } catch (...) {
    throw std::runtime_error("Topology: invalid double token for " + what + ": '" + x + "'");
  }
}
} // namespace detail_topo

// Read a simple bond table.
// Supported formats per non-comment line:
//   - "id_i id_j" (2 columns)
//   - "bond_type id_i id_j" (3 columns)
//   - "bond_id bond_type id_i id_j" (4+ columns; extra columns ignored; last two are atom ids)
inline Topology read_bond_table(const fs::path& path) {
  std::ifstream ifs(path);
  if (!ifs) throw std::runtime_error("Topology: failed to open bond table: " + path.string());

  Topology topo;
  topo.loaded_bonds = true;

  std::string line;
  while (std::getline(ifs, line)) {
    std::string s = detail_topo::trim(detail_topo::strip_comment(line));
    if (s.empty()) continue;
    const auto toks = detail_topo::split_ws(s);
    if (toks.size() < 2) continue;

    int btype = 0;
    std::int64_t ai = 0, aj = 0;
    if (toks.size() == 2) {
      ai = detail_topo::to_i64(toks[0], "bond_table.atom_i");
      aj = detail_topo::to_i64(toks[1], "bond_table.atom_j");
    } else if (toks.size() == 3) {
      btype = detail_topo::to_int(toks[0], "bond_table.bond_type");
      ai = detail_topo::to_i64(toks[1], "bond_table.atom_i");
      aj = detail_topo::to_i64(toks[2], "bond_table.atom_j");
    } else {
      // 4+ cols: assume bond_id bond_type ... atom_i atom_j
      btype = detail_topo::to_int(toks[1], "bond_table.bond_type");
      ai = detail_topo::to_i64(toks[toks.size() - 2], "bond_table.atom_i");
      aj = detail_topo::to_i64(toks[toks.size() - 1], "bond_table.atom_j");
    }
    topo.bond_ids.push_back(Topology::BondId{btype, ai, aj});
  }

  return topo;
}

// Read requested topology sections from a LAMMPS data file.
//
// Notes:
// - Sections are detected by header lines like "Masses", "Bonds", etc (optionally with trailing comments).
// - The reader is strict: if a requested section is not found, it throws.
inline Topology read_lammps_data_topology(const fs::path& path, const TopologySectionMask& mask) {
  std::ifstream ifs(path);
  if (!ifs) throw std::runtime_error("Topology: failed to open LAMMPS data file: " + path.string());

  enum class Sec { None, Masses, Bonds, Angles, Dihedrals, Impropers, Other };
  auto classify = [](const std::string& line) -> Sec {
    // We expect already-trimmed, comment-stripped input.
    if (line.rfind("Masses", 0) == 0) return Sec::Masses;
    if (line.rfind("Bonds", 0) == 0) return Sec::Bonds;
    if (line.rfind("Angles", 0) == 0) return Sec::Angles;
    if (line.rfind("Dihedrals", 0) == 0) return Sec::Dihedrals;
    if (line.rfind("Impropers", 0) == 0) return Sec::Impropers;
    return Sec::Other;
  };

  Topology topo;
  std::string line;
  Sec cur = Sec::None;

  bool found_masses = false;
  bool found_bonds = false;
  bool found_angles = false;
  bool found_dihedrals = false;
  bool found_impropers = false;

  while (std::getline(ifs, line)) {
    const std::string raw = detail_topo::trim(detail_topo::strip_comment(line));
    if (raw.empty()) continue;

    if (detail_topo::starts_with_alpha(raw)) {
      cur = classify(raw);
      if (cur == Sec::Masses) { found_masses = true; topo.loaded_masses = mask.masses; }
      if (cur == Sec::Bonds) { found_bonds = true; topo.loaded_bonds = mask.bonds; }
      if (cur == Sec::Angles) { found_angles = true; topo.loaded_angles = mask.angles; }
      if (cur == Sec::Dihedrals) { found_dihedrals = true; topo.loaded_dihedrals = mask.dihedrals; }
      if (cur == Sec::Impropers) { found_impropers = true; topo.loaded_impropers = mask.impropers; }
      continue;
    }

    // data line
    const auto toks = detail_topo::split_ws(raw);
    if (toks.empty()) continue;

    switch (cur) {
      case Sec::Masses: {
        if (!mask.masses) break;
        if (toks.size() < 2) throw std::runtime_error("Topology: malformed Masses line (need >=2 cols): " + raw);
        const int type = detail_topo::to_int(toks[0], "masses.type");
        const double m = detail_topo::to_f64(toks[1], "masses.mass");
        if (type <= 0) throw std::runtime_error("Topology: Masses type must be >= 1: " + raw);
        if (static_cast<std::size_t>(type) >= topo.mass_by_type.size()) {
          topo.mass_by_type.resize(static_cast<std::size_t>(type) + 1, 0.0);
        }
        topo.mass_by_type[static_cast<std::size_t>(type)] = m;
        break;
      }
      case Sec::Bonds: {
        if (!mask.bonds) break;
        if (toks.size() < 4) throw std::runtime_error("Topology: malformed Bonds line (need >=4 cols): " + raw);
        const int btype = detail_topo::to_int(toks[1], "bonds.bond_type");
        const std::int64_t ai = detail_topo::to_i64(toks[toks.size() - 2], "bonds.atom_i");
        const std::int64_t aj = detail_topo::to_i64(toks[toks.size() - 1], "bonds.atom_j");
        topo.bond_ids.push_back(Topology::BondId{btype, ai, aj});
        break;
      }
      case Sec::Angles: {
        if (!mask.angles) break;
        if (toks.size() < 5) throw std::runtime_error("Topology: malformed Angles line (need >=5 cols): " + raw);
        const int atype = detail_topo::to_int(toks[1], "angles.angle_type");
        const std::int64_t ai = detail_topo::to_i64(toks[toks.size() - 3], "angles.atom_i");
        const std::int64_t aj = detail_topo::to_i64(toks[toks.size() - 2], "angles.atom_j");
        const std::int64_t ak = detail_topo::to_i64(toks[toks.size() - 1], "angles.atom_k");
        topo.angle_ids.push_back(Topology::AngleId{atype, ai, aj, ak});
        break;
      }
      case Sec::Dihedrals: {
        if (!mask.dihedrals) break;
        if (toks.size() < 6) throw std::runtime_error("Topology: malformed Dihedrals line (need >=6 cols): " + raw);
        const int dtype = detail_topo::to_int(toks[1], "dihedrals.dihedral_type");
        const std::int64_t ai = detail_topo::to_i64(toks[toks.size() - 4], "dihedrals.atom_i");
        const std::int64_t aj = detail_topo::to_i64(toks[toks.size() - 3], "dihedrals.atom_j");
        const std::int64_t ak = detail_topo::to_i64(toks[toks.size() - 2], "dihedrals.atom_k");
        const std::int64_t al = detail_topo::to_i64(toks[toks.size() - 1], "dihedrals.atom_l");
        topo.dihedral_ids.push_back(Topology::DihedralId{dtype, ai, aj, ak, al});
        break;
      }
      case Sec::Impropers: {
        if (!mask.impropers) break;
        if (toks.size() < 6) throw std::runtime_error("Topology: malformed Impropers line (need >=6 cols): " + raw);
        const int itype = detail_topo::to_int(toks[1], "impropers.improper_type");
        const std::int64_t ai = detail_topo::to_i64(toks[toks.size() - 4], "impropers.atom_i");
        const std::int64_t aj = detail_topo::to_i64(toks[toks.size() - 3], "impropers.atom_j");
        const std::int64_t ak = detail_topo::to_i64(toks[toks.size() - 2], "impropers.atom_k");
        const std::int64_t al = detail_topo::to_i64(toks[toks.size() - 1], "impropers.atom_l");
        topo.improper_ids.push_back(Topology::ImproperId{itype, ai, aj, ak, al});
        break;
      }
      default:
        break;
    }
  }

  // Validate requested sections exist.
  if (mask.masses && !found_masses) throw std::runtime_error("Topology: LAMMPS data file has no 'Masses' section: " + path.string());
  if (mask.bonds && !found_bonds) throw std::runtime_error("Topology: LAMMPS data file has no 'Bonds' section: " + path.string());
  if (mask.angles && !found_angles) throw std::runtime_error("Topology: LAMMPS data file has no 'Angles' section: " + path.string());
  if (mask.dihedrals && !found_dihedrals) throw std::runtime_error("Topology: LAMMPS data file has no 'Dihedrals' section: " + path.string());
  if (mask.impropers && !found_impropers) throw std::runtime_error("Topology: LAMMPS data file has no 'Impropers' section: " + path.string());

  // For LAMMPS data, we mark sections loaded if the mask requested them and we found headers.
  // (Data lists may still be empty; that's allowed.)
  topo.loaded_masses = mask.masses && found_masses;
  topo.loaded_bonds = mask.bonds && found_bonds;
  topo.loaded_angles = mask.angles && found_angles;
  topo.loaded_dihedrals = mask.dihedrals && found_dihedrals;
  topo.loaded_impropers = mask.impropers && found_impropers;

  return topo;
}

} // namespace pilots
