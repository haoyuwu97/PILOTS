#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "pilots/alg/geom/GeometryView.hpp"
#include "pilots/alg/graph/Components.hpp"

namespace pilots::alg::graph {

// Non-PBC spanning check based on axis-aligned bounding boxes in *fractional* coordinates
// (lambda) inside the primary box.
//
// This is an intentionally conservative v1 criterion, useful for quick audits and
// non-periodic systems. For periodic spanning/percolation, a dedicated v1.1 module
// will be needed.
//
// - We wrap positions into the primary box before computing bounds.
// - Bounds are computed in lambda-space, which correctly represents the parallelepiped
//   faces for triclinic boxes.
// - A component is said to span an axis if its wrapped lambda-AABB touches both faces
//   (0 and 1) along that axis within a tolerance.
inline NonPBCSpanningResult compute_nonpbc_spanning(
    const pilots::Box& box,
    const pilots::alg::geom::GeometryView& geo,
    const ComponentsResult& comps,
    double eps_fraction = 1e-6) {

  if (geo.size() != comps.component_id.size()) {
    throw std::runtime_error("compute_nonpbc_spanning: geo.size != component_id.size");
  }

  const std::size_t C = comps.component_sizes.size();
  NonPBCSpanningResult out;
  out.spans_x.assign(C, false);
  out.spans_y.assign(C, false);
  out.spans_z.assign(C, false);

  if (C == 0) return out;

  const double eps = std::max(0.0, eps_fraction);

  std::vector<double> mins(C, 0.0), maxs(C, 0.0);
  std::vector<double> mint(C, 0.0), maxt(C, 0.0);
  std::vector<double> minu(C, 0.0), maxu(C, 0.0);
  std::vector<bool> init(C, false);

  for (std::size_t i = 0; i < geo.size(); ++i) {
    const std::size_t cid = comps.component_id[i];
    if (cid >= C) throw std::runtime_error("compute_nonpbc_spanning: bad component id");

    auto r = geo.wrap_position(geo.position(i));
    const auto lam = box.to_lambda(r[0], r[1], r[2]);
    const double s = lam[0], t = lam[1], u = lam[2];

    if (!init[cid]) {
      init[cid] = true;
      mins[cid] = maxs[cid] = s;
      mint[cid] = maxt[cid] = t;
      minu[cid] = maxu[cid] = u;
    } else {
      mins[cid] = std::min(mins[cid], s);
      maxs[cid] = std::max(maxs[cid], s);
      mint[cid] = std::min(mint[cid], t);
      maxt[cid] = std::max(maxt[cid], t);
      minu[cid] = std::min(minu[cid], u);
      maxu[cid] = std::max(maxu[cid], u);
    }
  }

  for (std::size_t cid = 0; cid < C; ++cid) {
    if (!init[cid]) continue;
    const bool sx = (mins[cid] <= eps) && (maxs[cid] >= 1.0 - eps);
    const bool sy = (mint[cid] <= eps) && (maxt[cid] >= 1.0 - eps);
    const bool sz = (minu[cid] <= eps) && (maxu[cid] >= 1.0 - eps);

    out.spans_x[cid] = sx;
    out.spans_y[cid] = sy;
    out.spans_z[cid] = sz;
    out.any_spans_x = out.any_spans_x || sx;
    out.any_spans_y = out.any_spans_y || sy;
    out.any_spans_z = out.any_spans_z || sz;
  }

  return out;
}

} // namespace pilots::alg::graph
