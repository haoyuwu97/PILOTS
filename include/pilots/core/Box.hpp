#pragma once

#include <array>
#include <cstdint>
#include <cmath>
#include <stdexcept>

namespace pilots {

struct Box {
  double xlo = 0.0, xhi = 0.0;
  double ylo = 0.0, yhi = 0.0;
  double zlo = 0.0, zhi = 0.0;
  // Tilt factors (triclinic)
  double xy = 0.0, xz = 0.0, yz = 0.0;
  bool triclinic = false;

  double lx() const { return xhi - xlo; }
  double ly() const { return yhi - ylo; }
  double lz() const { return zhi - zlo; }

  void validate() const {
    const double Lx = lx();
    const double Ly = ly();
    const double Lz = lz();
    if (Lx == 0.0 || Ly == 0.0 || Lz == 0.0) {
      throw std::runtime_error("Box: invalid box lengths (zero)");
    }
  }

  // Unwrap (x,y,z) using image flags. Works for orthorhombic and triclinic boxes.
  std::array<double,3> unwrap(double x, double y, double z,
                              std::int64_t ix, std::int64_t iy, std::int64_t iz) const {
    validate();

    const double Lx = lx();
    const double Ly = ly();
    const double Lz = lz();

    if (!triclinic) {
      return {x + static_cast<double>(ix) * Lx,
              y + static_cast<double>(iy) * Ly,
              z + static_cast<double>(iz) * Lz};
    }

    // LAMMPS triclinic unwrapping using tilt factors:
    // x' = x + ix*Lx + iy*xy + iz*xz
    // y' = y + iy*Ly + iz*yz
    // z' = z + iz*Lz
    return {x + static_cast<double>(ix) * Lx + static_cast<double>(iy) * xy + static_cast<double>(iz) * xz,
            y + static_cast<double>(iy) * Ly + static_cast<double>(iz) * yz,
            z + static_cast<double>(iz) * Lz};
  }

  // Convert a Cartesian point to fractional coordinates (lambda) in [0,1) without wrapping.
  // For triclinic: r = origin + H * lambda, H = [[Lx, xy, xz],[0, Ly, yz],[0,0,Lz]].
  std::array<double,3> to_lambda(double x, double y, double z) const {
    validate();

    const double rx = x - xlo;
    const double ry = y - ylo;
    const double rz = z - zlo;

    const double Lx = lx();
    const double Ly = ly();
    const double Lz = lz();

    if (!triclinic) {
      return {rx / Lx, ry / Ly, rz / Lz};
    }

    const double u = rz / Lz;
    const double t = (ry - u * yz) / Ly;
    const double s = (rx - t * xy - u * xz) / Lx;
    return {s, t, u};
  }

  // Convert fractional coordinates (lambda) to Cartesian.
  std::array<double,3> from_lambda(double s, double t, double u) const {
    validate();

    const double Lx = lx();
    const double Ly = ly();
    const double Lz = lz();

    if (!triclinic) {
      return {xlo + s * Lx, ylo + t * Ly, zlo + u * Lz};
    }

    const double x = xlo + s * Lx + t * xy + u * xz;
    const double y = ylo + t * Ly + u * yz;
    const double z = zlo + u * Lz;
    return {x, y, z};
  }

  // Wrap a Cartesian point back into the primary box (lambda in [0,1)).
  std::array<double,3> wrap(double x, double y, double z) const {
    auto lam = to_lambda(x, y, z);
    lam[0] -= std::floor(lam[0]);
    lam[1] -= std::floor(lam[1]);
    lam[2] -= std::floor(lam[2]);
    return from_lambda(lam[0], lam[1], lam[2]);
  }

  // Minimum-image displacement vector from r1 to r2 (r2 - r1) with PBC.
  std::array<double,3> min_image_displacement(double x1, double y1, double z1,
                                              double x2, double y2, double z2) const {
    validate();

    const double Lx = lx();
    const double Ly = ly();
    const double Lz = lz();

    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;

    if (!triclinic) {
      dx -= std::nearbyint(dx / Lx) * Lx;
      dy -= std::nearbyint(dy / Ly) * Ly;
      dz -= std::nearbyint(dz / Lz) * Lz;
      return {dx, dy, dz};
    }

    // Convert displacement to fractional components via inverse(H).
    const double u = dz / Lz;
    const double t = (dy - u * yz) / Ly;
    const double s = (dx - t * xy - u * xz) / Lx;

    double sw = s - std::nearbyint(s);
    double tw = t - std::nearbyint(t);
    double uw = u - std::nearbyint(u);

    dx = sw * Lx + tw * xy + uw * xz;
    dy = tw * Ly + uw * yz;
    dz = uw * Lz;

    return {dx, dy, dz};
  }

  std::array<double,3> min_image_displacement(const std::array<double,3>& r1,
                                              const std::array<double,3>& r2) const {
    return min_image_displacement(r1[0], r1[1], r1[2], r2[0], r2[1], r2[2]);
  }
};

} // namespace pilots
