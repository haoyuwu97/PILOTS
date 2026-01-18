#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pilots/core/Frame.hpp"
#include "pilots/select/CombineExpr.hpp"
#include "pilots/select/GroupRegistry.hpp"
#include "pilots/select/Selection.hpp"
#include "pilots/select/SelectionView.hpp"
#include "pilots/select/TopoGroupRegistry.hpp"
#include "pilots/topology/Topology.hpp"

namespace pilots {

// Shared selection resolver with caching.
//
// It unifies AtomGroup ([groups]) and TopoGroup ([topo_groups]) into atom sets,
// and supports measure-level combination via a boolean expression over variables A and T.
//
// C1/C2/C3 requirements in this package:
// - Two inputs (AtomGroup + TopoGroup), both producing atom sets
// - Combine expression over A,T with &,|,! and parentheses
// - Static selections cached once per unique spec
// - Dynamic selections cached per frame (by frame_index) and audited (size stats)
class SelectionProvider {
public:
  using SizeStats = GroupRegistry::SizeStats;

  SelectionProvider(GroupRegistry* atom_groups,
                    TopoGroupRegistry* topo_groups,
                    const Topology* topology)
  : atom_groups_(atom_groups), topo_groups_(topo_groups), topology_(topology) {
    if (!atom_groups_) throw std::runtime_error("SelectionProvider: atom_groups is null");
  }

  struct SelectionAuditInfo {
    std::string key;
    bool is_dynamic = false;
    std::size_t static_size = 0;
    SizeStats dynamic_size;
  };

  // Determine whether a combined selection spec is dynamic (depends on dynamic group/topo_group).
  bool is_dynamic_spec(const std::string& group_ref, const std::string& topo_group_ref) const {
    const bool dynA = is_dynamic_atom_ref_(group_ref);
    const bool dynT = is_dynamic_topo_ref_(topo_group_ref);
    return dynA || dynT;
  }

  // Resolve and return a cached view for the combined selection.
  //
  // If topo_group_ref is non-empty and not "all", topology must be present and bonds loaded.
  SelectionView get_combined_view(const Frame& frame,
                                 std::size_t frame_index,
                                 const std::string& group_ref,
                                 const std::string& topo_group_ref,
                                 const std::string& combine_expr) {
    const std::string g = group_ref.empty() ? std::string("all") : group_ref;
    const std::string t = topo_group_ref.empty() ? std::string("all") : topo_group_ref;
    const std::string key = make_key_(g, t, combine_expr);

    const bool dyn = is_dynamic_spec(g, t);
    if (!dyn) {
      auto it = static_combined_.find(key);
      if (it != static_combined_.end()) return view_of(it->second);
      // compute once and cache
      Selection out = compute_combined_(frame, frame_index, g, t, combine_expr);
      out.name = key;
      auto [ins, ok] = static_combined_.emplace(key, std::move(out));
      return view_of(ins->second);
    }

    auto& entry = dynamic_combined_[key];
    if (entry.last_frame_index == frame_index) {
      return view_of(entry.sel);
    }
    entry.sel = compute_combined_(frame, frame_index, g, t, combine_expr);
    entry.sel.name = key;
    entry.last_frame_index = frame_index;
    entry.size_stats.add(entry.sel.idx.size());
    return view_of(entry.sel);
  }

  std::vector<SelectionAuditInfo> audit() const {
    std::vector<SelectionAuditInfo> out;
    out.reserve(static_combined_.size() + dynamic_combined_.size());
    for (const auto& kv : static_combined_) {
      SelectionAuditInfo a;
      a.key = kv.first;
      a.is_dynamic = false;
      a.static_size = kv.second.idx.size();
      out.push_back(std::move(a));
    }
    for (const auto& kv : dynamic_combined_) {
      SelectionAuditInfo a;
      a.key = kv.first;
      a.is_dynamic = true;
      a.dynamic_size = kv.second.size_stats;
      out.push_back(std::move(a));
    }
    std::sort(out.begin(), out.end(), [](const SelectionAuditInfo& a, const SelectionAuditInfo& b) {
      return a.key < b.key;
    });
    return out;
  }

private:
  GroupRegistry* atom_groups_ = nullptr;
  TopoGroupRegistry* topo_groups_ = nullptr;
  const Topology* topology_ = nullptr;

  // Ad-hoc selector expressions (those not registered as named groups).
  // We treat them as *static* unless they contain a region clause, because regions
  // depend on the instantaneous box/positions.
  std::unordered_map<std::string, Selection> static_selector_cache_;

  std::unordered_map<std::string, Selection> static_combined_;

  struct DynSelEntry {
    Selection sel;
    std::size_t last_frame_index = static_cast<std::size_t>(-1);
  };
  // Dynamic ad-hoc selector expressions (currently: selectors that include a region).
  // Cached per-frame to avoid redundant recomputation within a frame.
  std::unordered_map<std::string, DynSelEntry> dynamic_selector_cache_;

  struct DynEntry {
    Selection sel;
    std::size_t last_frame_index = static_cast<std::size_t>(-1);
    SizeStats size_stats;
  };
  std::unordered_map<std::string, DynEntry> dynamic_combined_;

  static std::string canonicalize_combine_(const std::string& expr_raw) {
    // Canonicalize the combine expression for stable caching/audit keys.
    // - Remove all ASCII whitespace
    // - Normalize variable names to uppercase (a->A, t->T)
    std::string out;
    out.reserve(expr_raw.size());
    for (char c : expr_raw) {
      const unsigned char uc = static_cast<unsigned char>(c);
      if (uc == ' ' || uc == '\t' || uc == '\r' || uc == '\n') continue;
      if (c == 'a') c = 'A';
      if (c == 't') c = 'T';
      out.push_back(c);
    }
    if (out.empty()) out = "A&T";
    return out;
  }

  static std::string make_key_(const std::string& g, const std::string& t, const std::string& c) {
    // A stable key string that also serves as a human-auditable name.
    std::string key;
    key.reserve(g.size() + t.size() + c.size() + 64);
    key += "A="; key += g;
    key += ";T="; key += t;
    key += ";combine="; key += canonicalize_combine_(c);
    return key;
  }

  static bool selector_is_dynamic_(const std::string& selector_expr) {
    // Only selector expressions (those containing ':') can be parsed here.
    // The only selector feature that is intrinsically frame-dependent is 'region',
    // because it uses the instantaneous Box/positions (PBC + triclinic aware).
    const auto spec = detail::parse_selector_expression(selector_expr);
    return spec.region.kind != detail::Region::Kind::None;
  }

  bool is_dynamic_atom_ref_(const std::string& group_ref) const {
    if (group_ref.empty() || group_ref == "all") return false;
    if (atom_groups_->has(group_ref)) return atom_groups_->is_dynamic(group_ref);
    // Ad-hoc selector expression: dynamic only when it includes a region clause.
    if (group_ref.find(':') != std::string::npos) {
      return selector_is_dynamic_(group_ref);
    }
    // Otherwise, this is not a selector (likely a typo or unsupported inline boolean expression).
    // Let evaluation throw a focused error later.
    return false;
  }

  bool is_dynamic_topo_ref_(const std::string& topo_ref) const {
    if (topo_ref.empty() || topo_ref == "all") return false;
    if (!topo_groups_) return false;
    if (topo_groups_->has(topo_ref)) return topo_groups_->is_dynamic(topo_ref);
    // Unrecognized topo refs are handled as errors during evaluation.
    return false;
  }

  const Selection& resolve_atom_ref_(const Frame& frame, std::size_t frame_index, const std::string& group_ref) {
    if (group_ref.empty() || group_ref == "all") {
      return atom_groups_->get("all", frame, frame_index);
    }
    if (atom_groups_->has(group_ref)) {
      return atom_groups_->get(group_ref, frame, frame_index);
    }
    // Ad-hoc selector expression.
    if (group_ref.find(':') == std::string::npos) {
      throw std::runtime_error(
          "SelectionProvider: unknown atom group '" + group_ref + "' (not in [groups]) and not a selector expression");
    }

    const bool dyn = selector_is_dynamic_(group_ref);
    if (!dyn) {
      auto it = static_selector_cache_.find(group_ref);
      if (it != static_selector_cache_.end()) return it->second;
      Selection sel = build_selection_from_selector(frame, group_ref, group_ref);
      sel.name = group_ref;
      auto [ins, ok] = static_selector_cache_.emplace(group_ref, std::move(sel));
      (void)ok;
      return ins->second;
    }

    auto& entry = dynamic_selector_cache_[group_ref];
    if (entry.last_frame_index == frame_index) return entry.sel;
    entry.sel = build_selection_from_selector(frame, group_ref, group_ref);
    entry.sel.name = group_ref;
    entry.last_frame_index = frame_index;
    return entry.sel;
  }

  const Selection& resolve_topo_ref_(const Frame& frame,
                                    std::size_t frame_index,
                                    const std::string& topo_ref) {
    if (topo_ref.empty() || topo_ref == "all") {
      return atom_groups_->get("all", frame, frame_index);
    }
    if (!topology_) {
      throw std::runtime_error("SelectionProvider: topo_group requested but no topology loaded");
    }
    if (!topology_->has_bonds()) {
      throw std::runtime_error("SelectionProvider: topo_group requested but topology has no bonds loaded");
    }
    if (!topo_groups_) {
      throw std::runtime_error("SelectionProvider: topo_group requested but topo_groups registry is missing");
    }
    if (!topo_groups_->has(topo_ref)) {
      throw std::runtime_error("SelectionProvider: unknown topo_group: '" + topo_ref + "'");
    }
    return topo_groups_->get(topo_ref, frame, frame_index, *topology_, *atom_groups_);
  }

  Selection compute_combined_(const Frame& frame,
                             std::size_t frame_index,
                             const std::string& group_ref,
                             const std::string& topo_ref,
                             const std::string& combine_expr) {
    const Selection& A = resolve_atom_ref_(frame, frame_index, group_ref);
    const Selection& T = resolve_topo_ref_(frame, frame_index, topo_ref);
    const std::vector<std::size_t> idx = combine_sets(A.idx, T.idx, frame.natoms, combine_expr);
    Selection out;
    out.idx = idx;
    return out;
  }
};

} // namespace pilots
