#pragma once

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "pilots/core/Frame.hpp"
#include "pilots/select/GroupRegistry.hpp"
#include "pilots/select/Selection.hpp"
#include "pilots/select/SelectionView.hpp"
#include "pilots/topology/Topology.hpp"

namespace pilots {

// Topology-backed atom groups.
//
// A TopoGroup always outputs an atom index set (Selection), but its definition
// can depend on loaded topology (bonds/graph) and/or an AtomGroup.
//
// This is intentionally NOT a general derived-property DSL: it is a small set
// of topology -> atom-set selectors.
class TopoGroupRegistry {
public:
  using SizeStats = GroupRegistry::SizeStats;

  struct TopoGroupDef {
    enum class Kind {
      AtomsInAnyBond,
      AtomsInBondType,
      AtomsInComponentSizeGe,
      AtomsIncidentToAtomGroup,
      AtomsDegreeIs,
      AtomsIsEnd,
    } kind = Kind::AtomsInAnyBond;

    bool dynamic = false;
    std::string expr; // raw expression (post selection_mode stripping)

    // Parameters for each kind.
    std::vector<int> bond_types;
    std::size_t component_size_ge = 0;
    std::string atom_group_ref; // for incident_to_group
    int degree_k = -1;
  };

  explicit TopoGroupRegistry(std::unordered_map<std::string, std::string> raw_defs)
  : raw_defs_(std::move(raw_defs)) {
    parse_all_();
  }

  bool has(const std::string& name) const {
    return defs_.find(name) != defs_.end();
  }

  bool is_dynamic(const std::string& name) const {
    auto it = defs_.find(name);
    if (it == defs_.end()) return false;
    return it->second.dynamic;
  }

  void prepare_static(const Frame& frame0, const Topology& topo, GroupRegistry& groups) {
    prepared_ = true;

    // built-in 'all'
    static_cache_["all"] = make_all_selection(frame0, "all");

    // Evaluate all non-dynamic topo groups once.
    for (const auto& kv : defs_) {
      const std::string& name = kv.first;
      if (name == "all") continue;
      const TopoGroupDef& def = kv.second;
      if (!def.dynamic) {
        // A static topo_group may refer to an AtomGroup (e.g. atoms_incident_to_group).
        // If the referenced AtomGroup is dynamic, using it here would silently freeze it
        // at frame 0, which is scientifically incorrect. Fail-fast.
        if (def.kind == TopoGroupDef::Kind::AtomsIncidentToAtomGroup && !def.atom_group_ref.empty()) {
          if (groups.has(def.atom_group_ref) && groups.is_dynamic(def.atom_group_ref)) {
            throw std::runtime_error(
                "TopoGroupRegistry: topo_group '" + name + "' is static but depends on dynamic atom group '" + def.atom_group_ref + "'. "
                "Add 'selection_mode:dynamic' to the topo_group, or make the atom group static.");
          }
        }
        static_cache_[name] = eval_(name, def, frame0, 0, topo, groups);
      }
    }
  }

  const Selection& get(const std::string& name,
                       const Frame& frame,
                       std::size_t frame_index,
                       const Topology& topo,
                       GroupRegistry& groups) {
    if (!prepared_) {
      throw std::runtime_error("TopoGroupRegistry: prepare_static() must be called before get()");
    }

    if (name == "all") {
      auto it = static_cache_.find("all");
      if (it == static_cache_.end()) throw std::runtime_error("TopoGroupRegistry: missing built-in 'all'");
      return it->second;
    }

    auto it = defs_.find(name);
    if (it == defs_.end()) {
      throw std::runtime_error("TopoGroupRegistry: unknown topo_group name: '" + name + "'");
    }
    const TopoGroupDef& def = it->second;
    if (!def.dynamic) {
      auto sit = static_cache_.find(name);
      if (sit == static_cache_.end()) {
        throw std::runtime_error("TopoGroupRegistry: static topo_group not prepared: '" + name + "'");
      }
      return sit->second;
    }

    auto& entry = dynamic_cache_[name];
    if (entry.last_frame_index == frame_index) {
      return entry.sel;
    }
    entry.sel = eval_(name, def, frame, frame_index, topo, groups);
    entry.last_frame_index = frame_index;
    entry.size_stats.add(entry.sel.idx.size());
    return entry.sel;
  }

  SelectionView get_view(const std::string& name,
                         const Frame& frame,
                         std::size_t frame_index,
                         const Topology& topo,
                         GroupRegistry& groups) {
    const Selection& s = get(name, frame, frame_index, topo, groups);
    return view_of(s);
  }

  struct TopoGroupAuditInfo {
    std::string name;
    std::string expr;
    bool is_dynamic = false;
    std::size_t static_size = 0;
    SizeStats size_stats; // samples=0 => never evaluated
  };

  std::vector<TopoGroupAuditInfo> audit() const {
    std::vector<TopoGroupAuditInfo> out;
    out.reserve(defs_.size() + 1);

    // built-in 'all'
    {
      TopoGroupAuditInfo g;
      g.name = "all";
      g.expr = "all";
      g.is_dynamic = false;
      if (auto it = static_cache_.find("all"); it != static_cache_.end()) {
        g.static_size = it->second.idx.size();
      }
      out.push_back(std::move(g));
    }

    for (const auto& kv : defs_) {
      const std::string& name = kv.first;
      if (name == "all") continue;
      const TopoGroupDef& def = kv.second;
      TopoGroupAuditInfo g;
      g.name = name;
      g.expr = def.expr;
      g.is_dynamic = def.dynamic;
      if (!def.dynamic) {
        if (auto it = static_cache_.find(name); it != static_cache_.end()) {
          g.static_size = it->second.idx.size();
        }
      } else {
        if (auto it = dynamic_cache_.find(name); it != dynamic_cache_.end()) {
          g.size_stats = it->second.size_stats;
        }
      }
      out.push_back(std::move(g));
    }

    std::sort(out.begin(), out.end(), [](const TopoGroupAuditInfo& a, const TopoGroupAuditInfo& b) {
      return a.name < b.name;
    });
    return out;
  }

private:
  struct DynEntry {
    Selection sel;
    std::size_t last_frame_index = static_cast<std::size_t>(-1);
    SizeStats size_stats;
  };

  std::unordered_map<std::string, std::string> raw_defs_;
  std::unordered_map<std::string, TopoGroupDef> defs_;

  bool prepared_ = false;
  std::unordered_map<std::string, Selection> static_cache_;
  std::unordered_map<std::string, DynEntry> dynamic_cache_;

  static std::string trim_(const std::string& s) {
    std::size_t b = 0;
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    std::size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
  }

  static std::vector<std::string> split_semicolon_(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char ch : s) {
      if (ch == ';') {
        auto t = trim_(cur);
        if (!t.empty()) out.push_back(t);
        cur.clear();
      } else {
        cur.push_back(ch);
      }
    }
    auto t = trim_(cur);
    if (!t.empty()) out.push_back(t);
    return out;
  }

  static std::string join_semicolon_(const std::vector<std::string>& parts, std::size_t from) {
    std::string out;
    for (std::size_t i = from; i < parts.size(); ++i) {
      if (!out.empty()) out += "; ";
      out += parts[i];
    }
    return out;
  }

  static std::vector<int> parse_int_list_(std::string s) {
    // Accept comma-separated with optional whitespace and/or braces.
    s = trim_(s);
    if (!s.empty() && s.front() == '{' && s.back() == '}') {
      s = s.substr(1, s.size() - 2);
    }
    std::vector<int> out;
    std::string cur;
    auto flush = [&]() {
      auto t = trim_(cur);
      if (!t.empty()) {
        std::size_t pos = 0;
        long long v = 0;
        try {
          v = std::stoll(t, &pos);
          if (pos != t.size()) throw std::invalid_argument("trailing");
        } catch (...) {
          throw std::runtime_error("TopoGroupRegistry: invalid integer in list: '" + t + "'");
        }
        out.push_back(static_cast<int>(v));
      }
      cur.clear();
    };
    for (char ch : s) {
      if (ch == ',' || ch == ' ' || ch == '\t') {
        flush();
      } else {
        cur.push_back(ch);
      }
    }
    flush();
    return out;
  }

  static std::string to_lower_(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
  }

  void parse_all_() {
    defs_.clear();
    // User defs
    for (const auto& kv : raw_defs_) {
      const std::string& name = kv.first;
      const std::string& value = kv.second;
      if (name.empty()) continue;
      if (name == "all") {
        // Disallow overriding built-in 'all' in topo_groups for now; it is always full selection.
        continue;
      }
      TopoGroupDef def;
      parse_one_(def, value);
      defs_[name] = std::move(def);
    }
  }

  void parse_one_(TopoGroupDef& def, const std::string& raw) {
    auto parts = split_semicolon_(raw);
    def.dynamic = false;

    std::size_t i = 0;
    for (; i < parts.size(); ++i) {
      const std::string p = trim_(parts[i]);
      if (p.rfind("selection_mode:", 0) == 0) {
        std::string mode = trim_(p.substr(std::string("selection_mode:").size()));
        mode = to_lower_(mode);
        if (mode == "dynamic") def.dynamic = true;
        else if (mode == "static") def.dynamic = false;
        else throw std::runtime_error("TopoGroupRegistry: invalid selection_mode: '" + mode + "'");
      } else {
        break;
      }
    }
    def.expr = join_semicolon_(parts, i);
    if (def.expr.empty()) throw std::runtime_error("TopoGroupRegistry: empty topo_group expression");

    // Parse kind + params.
    const std::string expr_l = to_lower_(trim_(def.expr));

    if (expr_l == "atoms_in_any_bond") {
      def.kind = TopoGroupDef::Kind::AtomsInAnyBond;
      return;
    }
    if (expr_l.rfind("atoms_in_bond_type", 0) == 0) {
      def.kind = TopoGroupDef::Kind::AtomsInBondType;
      std::string rest = trim_(def.expr.substr(std::string("atoms_in_bond_type").size()));
      if (!rest.empty() && rest.front() == ':') rest = trim_(rest.substr(1));
      def.bond_types = parse_int_list_(rest);
      if (def.bond_types.empty()) throw std::runtime_error("TopoGroupRegistry: atoms_in_bond_type requires a non-empty list");
      return;
    }
    if (expr_l.rfind("atoms_in_component_size", 0) == 0) {
      // Accept: atoms_in_component_size >= N  OR atoms_in_component_size>=N
      std::string rest = trim_(def.expr.substr(std::string("atoms_in_component_size").size()));
      rest = trim_(rest);
      // remove spaces
      std::string rs;
      for (char c : rest) if (!std::isspace(static_cast<unsigned char>(c))) rs.push_back(c);
      if (rs.rfind(">=", 0) != 0) {
        throw std::runtime_error("TopoGroupRegistry: atoms_in_component_size requires '>= N'" );
      }
      std::string n_s = rs.substr(2);
      if (n_s.empty()) throw std::runtime_error("TopoGroupRegistry: atoms_in_component_size missing N");
      std::size_t pos = 0;
      unsigned long long v = 0;
      try {
        v = std::stoull(n_s, &pos);
        if (pos != n_s.size()) throw std::invalid_argument("trailing");
      } catch (...) {
        throw std::runtime_error("TopoGroupRegistry: atoms_in_component_size invalid N: '" + n_s + "'");
      }
      def.kind = TopoGroupDef::Kind::AtomsInComponentSizeGe;
      def.component_size_ge = static_cast<std::size_t>(v);
      return;
    }
    if (expr_l.rfind("atoms_incident_to_group", 0) == 0 || expr_l.rfind("atoms_incident_to_atom_group", 0) == 0) {
      def.kind = TopoGroupDef::Kind::AtomsIncidentToAtomGroup;
      // Allow atoms_incident_to_group(G) or atoms_incident_to_group:G
      std::string rest = def.expr;
      // find first '(' or ':' after the keyword
      const std::string k1 = "atoms_incident_to_group";
      const std::string k2 = "atoms_incident_to_atom_group";
      std::size_t off = std::string::npos;
      if (expr_l.rfind(k2,0)==0) off = k2.size();
      else off = k1.size();
      rest = trim_(def.expr.substr(off));
      if (!rest.empty() && rest.front() == ':') rest = trim_(rest.substr(1));
      if (!rest.empty() && rest.front() == '(' && rest.back() == ')') {
        rest = trim_(rest.substr(1, rest.size()-2));
      }
      if (rest.empty()) throw std::runtime_error("TopoGroupRegistry: atoms_incident_to_group requires a group name");
      def.atom_group_ref = rest;
      // This topo group is dynamic if the referenced atom group is dynamic.
      // We cannot resolve this here without GroupRegistry; Runner should validate.
      return;
    }
    if (expr_l.rfind("atoms_degree_is", 0) == 0) {
      def.kind = TopoGroupDef::Kind::AtomsDegreeIs;
      std::string rest = trim_(def.expr.substr(std::string("atoms_degree_is").size()));
      if (!rest.empty() && rest.front() == ':') rest = trim_(rest.substr(1));
      if (rest.empty()) throw std::runtime_error("TopoGroupRegistry: atoms_degree_is missing k");
      std::size_t pos = 0;
      long long v = 0;
      try {
        v = std::stoll(rest, &pos);
        if (pos != rest.size()) throw std::invalid_argument("trailing");
      } catch (...) {
        throw std::runtime_error("TopoGroupRegistry: atoms_degree_is invalid k: '" + rest + "'");
      }
      if (v < 0) throw std::runtime_error("TopoGroupRegistry: atoms_degree_is requires k >= 0");
      def.degree_k = static_cast<int>(v);
      return;
    }
    if (expr_l == "atoms_is_end") {
      def.kind = TopoGroupDef::Kind::AtomsIsEnd;
      return;
    }

    throw std::runtime_error("TopoGroupRegistry: unknown topo_group selector: '" + def.expr + "'");
  }

  Selection eval_(const std::string& name,
                 const TopoGroupDef& def,
                 const Frame& frame,
                 std::size_t frame_index,
                 const Topology& topo,
                 GroupRegistry& groups) {
    if (!topo.finalized) {
      throw std::runtime_error("TopoGroupRegistry: topology must be finalized before evaluating topo_groups");
    }
    if (!topo.has_bonds() && def.kind != TopoGroupDef::Kind::AtomsInAnyBond) {
      // All current topo selectors require bonds/graph except the 'all' group.
      // atoms_in_any_bond also needs bonds, but we check separately for better messaging.
    }
    if (!topo.has_bonds()) {
      // atoms_in_any_bond and everything else requires bonds.
      throw std::runtime_error("TopoGroupRegistry: topo_groups require bonds; no bonds loaded in topology");
    }

    std::vector<unsigned char> mark(topo.natoms, 0);

    switch (def.kind) {
      case TopoGroupDef::Kind::AtomsInAnyBond: {
        for (const auto& b : topo.bonds) {
          mark[b.i] = 1;
          mark[b.j] = 1;
        }
        break;
      }
      case TopoGroupDef::Kind::AtomsInBondType: {
        // Bond type ids can be sparse (large max with few used types). Use a dense
        // table only when it is safe; otherwise fall back to a set.
        int tmax = -1;
        for (int t : def.bond_types) if (t > tmax) tmax = t;
        const std::size_t ntypes = def.bond_types.size();
        const bool use_dense = (tmax >= 0) &&
                               (static_cast<std::size_t>(tmax) <= 200000u) &&
                               (static_cast<std::size_t>(tmax) <= ntypes * 16u);

        if (use_dense) {
          std::vector<unsigned char> allow(static_cast<std::size_t>(tmax + 1), 0);
          for (int t : def.bond_types) {
            if (t >= 0 && static_cast<std::size_t>(t) < allow.size()) allow[static_cast<std::size_t>(t)] = 1;
          }
          for (const auto& b : topo.bonds) {
            const int t = b.type;
            if (t >= 0 && static_cast<std::size_t>(t) < allow.size() && allow[static_cast<std::size_t>(t)]) {
              mark[b.i] = 1;
              mark[b.j] = 1;
            }
          }
        } else {
          std::unordered_set<int> allow(def.bond_types.begin(), def.bond_types.end());
          for (const auto& b : topo.bonds) {
            if (allow.find(b.type) != allow.end()) {
              mark[b.i] = 1;
              mark[b.j] = 1;
            }
          }
        }
        break;
      }
      case TopoGroupDef::Kind::AtomsInComponentSizeGe: {
        topo.ensure_components();
        for (std::size_t i = 0; i < topo.natoms; ++i) {
          const std::size_t cid = topo.component_id[i];
          if (cid < topo.component_sizes.size() && topo.component_sizes[cid] >= def.component_size_ge) {
            mark[i] = 1;
          }
        }
        break;
      }
      case TopoGroupDef::Kind::AtomsIncidentToAtomGroup: {
        const Selection& gsel = groups.get(def.atom_group_ref, frame, frame_index);
        for (const std::size_t u : gsel.idx) {
          if (u >= topo.adjacency.size()) continue;
          for (const std::size_t v : topo.adjacency[u]) {
            mark[v] = 1;
          }
        }
        break;
      }
      case TopoGroupDef::Kind::AtomsDegreeIs: {
        const int k = def.degree_k;
        for (std::size_t i = 0; i < topo.natoms; ++i) {
          const std::size_t deg = (i < topo.adjacency.size()) ? topo.adjacency[i].size() : 0;
          if (static_cast<int>(deg) == k) mark[i] = 1;
        }
        break;
      }
      case TopoGroupDef::Kind::AtomsIsEnd: {
        for (std::size_t i = 0; i < topo.natoms; ++i) {
          const std::size_t deg = (i < topo.adjacency.size()) ? topo.adjacency[i].size() : 0;
          if (deg == 1) mark[i] = 1;
        }
        break;
      }
    }

    Selection out;
    out.name = name;
    out.idx.reserve(topo.natoms / 4);
    for (std::size_t i = 0; i < topo.natoms; ++i) {
      if (mark[i]) out.idx.push_back(i);
    }
    return out;
  }
};

} // namespace pilots
