#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pilots/core/Frame.hpp"

namespace pilots {

struct Selection {
  std::string name;
  std::vector<std::size_t> idx; // indices into canonical atom arrays
};

inline Selection make_all_selection(const Frame& frame, std::string name = "all") {
  Selection sel;
  sel.name = std::move(name);
  sel.idx.resize(frame.natoms);
  for (std::size_t i = 0; i < frame.natoms; ++i) sel.idx[i] = i;
  return sel;
}

namespace detail {

inline std::string trim(const std::string& in) {
  std::size_t b = 0;
  while (b < in.size() && (in[b] == ' ' || in[b] == '\t' || in[b] == '\r' || in[b] == '\n')) ++b;
  std::size_t e = in.size();
  while (e > b && (in[e - 1] == ' ' || in[e - 1] == '\t' || in[e - 1] == '\r' || in[e - 1] == '\n')) --e;
  return in.substr(b, e - b);
}

inline std::vector<std::string> split_by_char(const std::string& s, char sep) {
  std::vector<std::string> out;
  std::string cur;
  for (char ch : s) {
    if (ch == sep) {
      auto t = trim(cur);
      if (!t.empty()) out.push_back(t);
      cur.clear();
    } else {
      cur.push_back(ch);
    }
  }
  auto t = trim(cur);
  if (!t.empty()) out.push_back(t);
  return out;
}

inline std::vector<std::string> split_csv_trim(const std::string& s) {
  return split_by_char(s, ',');
}

inline void parse_int_or_range(const std::string& token,
                              std::vector<std::pair<std::int64_t, std::int64_t>>& ranges) {
  // token: "a" or "a-b" (inclusive)
  auto dash = token.find('-');
  if (dash == std::string::npos) {
    std::int64_t v = 0;
    try {
      std::size_t pos = 0;
      v = std::stoll(token, &pos);
      if (pos != token.size()) throw std::invalid_argument("trailing");
    } catch (...) {
      throw std::runtime_error("Selection: invalid integer token: '" + token + "'");
    }
    ranges.emplace_back(v, v);
    return;
  }

  const std::string a = token.substr(0, dash);
  const std::string b = token.substr(dash + 1);
  if (a.empty() || b.empty()) {
    throw std::runtime_error("Selection: invalid range token: '" + token + "'");
  }
  std::int64_t va = 0, vb = 0;
  try {
    std::size_t pos1 = 0, pos2 = 0;
    va = std::stoll(a, &pos1);
    vb = std::stoll(b, &pos2);
    if (pos1 != a.size() || pos2 != b.size()) throw std::invalid_argument("trailing");
  } catch (...) {
    throw std::runtime_error("Selection: invalid range token: '" + token + "'");
  }
  if (vb < va) std::swap(va, vb);
  ranges.emplace_back(va, vb);
}

inline bool in_ranges(std::int64_t v, const std::vector<std::pair<std::int64_t,std::int64_t>>& ranges) {
  for (auto [a,b] : ranges) {
    if (v >= a && v <= b) return true;
  }
  return false;
}

struct Region {
  enum class Kind { None, Box, Slab, Sphere } kind = Kind::None;
  // Box: [xmin,xmax]x[ymin,ymax]x[zmin,zmax]
  double xmin=0, xmax=0, ymin=0, ymax=0, zmin=0, zmax=0;
  // Slab: axis in {x,y,z} with [min,max]
  char axis='z';
  double amin=0, amax=0;
  // Sphere: center and radius
  double cx=0, cy=0, cz=0, r=0;

  bool contains(const Box& box, double x, double y, double z) const {
    switch (kind) {
      case Kind::None: return true;
      case Kind::Box: {
        const auto w = box.wrap(x, y, z);
        const double wx = w[0], wy = w[1], wz = w[2];
        return (wx >= xmin && wx <= xmax && wy >= ymin && wy <= ymax && wz >= zmin && wz <= zmax);
      }
      case Kind::Slab: {
        const auto w = box.wrap(x, y, z);
        const double v = (axis=='x')?w[0]:((axis=='y')?w[1]:w[2]);
        return (v >= amin && v <= amax);
      }
      case Kind::Sphere: {
        const auto cw = box.wrap(cx, cy, cz);
        const auto pw = box.wrap(x, y, z);
        const auto d = box.min_image_displacement(cw, pw);
        return (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]) <= r*r;
      }
    }
    return true;
  }
};

struct SelectorSpec {
  bool all = false;
  bool has_type = false;
  bool has_mol  = false;
  bool has_id   = false;
  std::vector<std::pair<std::int64_t,std::int64_t>> type_ranges;
  std::vector<std::pair<std::int64_t,std::int64_t>> mol_ranges;
  std::vector<std::pair<std::int64_t,std::int64_t>> id_ranges;
  Region region;
};

inline double parse_double_strict(const std::string& s, const std::string& what) {
  try {
    std::size_t pos = 0;
    double v = std::stod(s, &pos);
    if (pos != s.size()) throw std::invalid_argument("trailing");
    return v;
  } catch (...) {
    throw std::runtime_error("Selection: invalid " + what + " value: '" + s + "'");
  }
}

inline SelectorSpec parse_selector_expression(const std::string& selector_raw) {
  SelectorSpec spec;
  std::string s = trim(selector_raw);
  if (s.empty() || s == "*" || s == "all") {
    spec.all = true;
    return spec;
  }

  // Allow compound clauses separated by ';'
  const auto clauses = split_by_char(s, ';');
  for (const auto& clause_raw : clauses) {
    const std::string clause = trim(clause_raw);
    if (clause.empty()) continue;

    auto colon = clause.find(':');
    if (colon == std::string::npos) {
      throw std::runtime_error("Selection: invalid clause (missing ':'): '" + clause + "'");
    }
    const std::string field = trim(clause.substr(0, colon));
    const std::string rhs   = trim(clause.substr(colon + 1));

    if (rhs == "*" || rhs == "all") {
      // means "no filter" for this field
      continue;
    }

    if (field == "type") {
      spec.has_type = true;
      for (const auto& item : split_csv_trim(rhs)) {
        parse_int_or_range(item, spec.type_ranges);
      }
    } else if (field == "mol" || field == "molecule") {
      spec.has_mol = true;
      for (const auto& item : split_csv_trim(rhs)) {
        parse_int_or_range(item, spec.mol_ranges);
      }
    } else if (field == "id") {
      spec.has_id = true;
      for (const auto& item : split_csv_trim(rhs)) {
        parse_int_or_range(item, spec.id_ranges);
      }
    } else if (field == "region") {
      // region:box,xmin,xmax,ymin,ymax,zmin,zmax
      // region:slab,axis,min,max
      // region:sphere,cx,cy,cz,r
      const auto parts = split_csv_trim(rhs);
      if (parts.empty()) throw std::runtime_error("Selection: empty region spec");
      const std::string kind = parts[0];
      Region r;
      if (kind == "box") {
        if (parts.size() != 7) {
          throw std::runtime_error("Selection: region:box expects 6 numbers: box,xmin,xmax,ymin,ymax,zmin,zmax");
        }
        r.kind = Region::Kind::Box;
        r.xmin = parse_double_strict(parts[1], "xmin");
        r.xmax = parse_double_strict(parts[2], "xmax");
        r.ymin = parse_double_strict(parts[3], "ymin");
        r.ymax = parse_double_strict(parts[4], "ymax");
        r.zmin = parse_double_strict(parts[5], "zmin");
        r.zmax = parse_double_strict(parts[6], "zmax");
        if (r.xmax < r.xmin) std::swap(r.xmin, r.xmax);
        if (r.ymax < r.ymin) std::swap(r.ymin, r.ymax);
        if (r.zmax < r.zmin) std::swap(r.zmin, r.zmax);
      } else if (kind == "slab") {
        if (parts.size() != 4) {
          throw std::runtime_error("Selection: region:slab expects 3 values: slab,axis,min,max");
        }
        r.kind = Region::Kind::Slab;
        if (parts[1].size() != 1 || (parts[1][0] != 'x' && parts[1][0] != 'y' && parts[1][0] != 'z')) {
          throw std::runtime_error("Selection: region:slab axis must be x,y,z");
        }
        r.axis = parts[1][0];
        r.amin = parse_double_strict(parts[2], "min");
        r.amax = parse_double_strict(parts[3], "max");
        if (r.amax < r.amin) std::swap(r.amin, r.amax);
      } else if (kind == "sphere") {
        if (parts.size() != 5) {
          throw std::runtime_error("Selection: region:sphere expects 5 values: sphere,cx,cy,cz,r");
        }
        r.kind = Region::Kind::Sphere;
        r.cx = parse_double_strict(parts[1], "cx");
        r.cy = parse_double_strict(parts[2], "cy");
        r.cz = parse_double_strict(parts[3], "cz");
        r.r  = parse_double_strict(parts[4], "r");
        if (r.r < 0) throw std::runtime_error("Selection: region:sphere r must be >= 0");
      } else {
        throw std::runtime_error("Selection: unsupported region kind: '" + kind + "'");
      }
      spec.region = r;
    } else {
      throw std::runtime_error("Selection: unsupported selector field: '" + field + "'");
    }
  }

  return spec;
}

} // namespace detail

// Selector grammar (extendable, config-friendly):
//
// Basic forms (backward compatible):
//   "*" or "all"                  -> all atoms
//   "type:*" or "type:1,2"        -> by type
//   "mol:1-100"                    -> by molecule id
//   "id:5,9-12"                    -> by atom id
//
// Compound clauses (AND semantics) separated by ';':
//   "type:1; region:slab,z,0,50"   -> type filter AND region filter
//
// Region clause:
//   region:box,xmin,xmax,ymin,ymax,zmin,zmax
//   region:slab,axis,min,max
//   region:sphere,cx,cy,cz,r
inline Selection build_selection_from_selector(const Frame& frame, std::string name, const std::string& selector) {
  Selection sel;
  sel.name = std::move(name);

  const auto spec = detail::parse_selector_expression(selector);
  if (spec.all) {
    return make_all_selection(frame, sel.name);
  }

  // Strict dependency checks for core metadata.
  // We only require type/mol if the selector uses them.
  const auto ids_all = frame.require_i64field("id");
  std::span<const int> types;
  std::span<const std::int64_t> mols;
  std::span<const std::int64_t> ids;
  if (spec.has_type) types = frame.require_intfield("type");
  if (spec.has_mol)  mols  = frame.require_i64field("mol");
  if (spec.has_id)   ids   = frame.require_i64field("id");

  // Coordinates are required for region filtering.
  const auto xu = frame.require_dfield("xu");
  const auto yu = frame.require_dfield("yu");
  const auto zu = frame.require_dfield("zu");

  sel.idx.reserve(frame.natoms);
  for (std::size_t i = 0; i < frame.natoms; ++i) {
    const std::int64_t idv = ids_all[i];
    const std::int64_t molv = spec.has_mol ? mols[i] : 0;
    const std::int64_t typev = spec.has_type ? static_cast<std::int64_t>(types[i]) : 0;

    if (spec.has_type && !detail::in_ranges(typev, spec.type_ranges)) continue;
    if (spec.has_mol  && !detail::in_ranges(molv,  spec.mol_ranges))  continue;
    if (spec.has_id   && !detail::in_ranges(idv,   spec.id_ranges))   continue;

    const double x = xu[i];
    const double y = yu[i];
    const double z = zu[i];
    if (!spec.region.contains(frame.box, x, y, z)) continue;

    sel.idx.push_back(i);
  }

  return sel;
}

// Build named groups from a map {name -> selector}. Always inject an "all" group.
inline std::unordered_map<std::string, Selection>
build_group_map(const Frame& frame, const std::unordered_map<std::string, std::string>& definitions) {
  std::unordered_map<std::string, Selection> groups;
  groups.reserve(definitions.size() + 1);
  groups.emplace("all", make_all_selection(frame, "all"));

  for (const auto& [name, sel] : definitions) {
    if (name == "all") {
      groups["all"] = build_selection_from_selector(frame, "all", sel);
      continue;
    }
    groups.emplace(name, build_selection_from_selector(frame, name, sel));
  }

  return groups;
}

} // namespace pilots
