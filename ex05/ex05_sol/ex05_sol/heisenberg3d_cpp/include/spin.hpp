#pragma once
// Structure to store a 3D spin.

#include <random>

#include "typedefs.hpp"

class Spin {
public:
  Spin() = default;
  Spin(Real x, Real y, Real z) : x_(x), y_(y), z_(z) {}
  Spin(const Spin& rhs) = default;
  Spin& operator=(const Spin& rhs) = default;
  Spin(Spin&& rhs) = default;
  Spin& operator=(Spin&& rhs) = default;

  template <class Rng>
  static Spin random(Rng& rng) {
    // std::uniform_real_distribution<Real> distro(-1., 1.); // This does not sample uniformly over a 3-sphere!
    std::normal_distribution<Real> distro(0., 1.); // use normal distribution to sample uniformly over 3-sphere
    Spin s(distro(rng), distro(rng), distro(rng));

    const Real norm_val = std::sqrt(s.norm2());
    s.x_ /= norm_val;
    s.y_ /= norm_val;
    s.z_ /= norm_val;
    return s;
  }

  Spin operator+(const Spin& rhs) const {
    return Spin(x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_);
  }

  Spin operator-(const Spin& rhs) const {
    return Spin(x_ - rhs.x_, y_ - rhs.y_, z_ - rhs.z_);
  }

  Spin operator-() const {
    return -1. * (*this);
  }

  Real operator*(const Spin& rhs) {
    return x_ * rhs.x_ + y_ * rhs.y_ + z_ * rhs.z_;
  }

  Real norm2() const {
    return x_ * x_ + y_ * y_ + z_ * z_;
  }

  friend inline Spin operator*(Real a, const Spin& s);

  // reflect the spin across the plane defined by random normal vector r:
  void flip(const Spin& r) { 
    const Spin delta = (2. * (*this * r)) * r;
    *this = *this - delta;
  }

private:
  Real x_;
  Real y_;
  Real z_;
  // TODO consider padding for vectorization.
};

Spin operator*(Real a, const Spin& s) {
  return Spin(a * s.x_, a * s.y_, a * s.z_);
}
