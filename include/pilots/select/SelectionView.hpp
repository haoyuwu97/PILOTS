#pragma once

#include <cstddef>
#include <span>
#include <string_view>

#include "pilots/select/Selection.hpp"

namespace pilots {

// A lightweight view over a Selection's index list.
//
// Lifetime: the underlying idx data must outlive the view.
// - For groups/topo_groups, this is typically owned by the (Topo)GroupRegistry caches.
// - For combined selections, this is owned by SelectionProvider caches.
struct SelectionView {
  std::string_view name;
  std::span<const std::size_t> idx;

  std::size_t size() const { return idx.size(); }
  bool empty() const { return idx.empty(); }
};

inline SelectionView view_of(const Selection& s) {
  return SelectionView{std::string_view(s.name), std::span<const std::size_t>(s.idx.data(), s.idx.size())};
}

} // namespace pilots
