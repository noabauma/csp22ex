/* Computational Statistical Physics 2021
 * Exercise 8
 * Author: Pascal Engeler, engelerp@phys.ethz.ch
 */
#ifndef LENNARD_JONES_HPP_INCLUDED
#define LENNARD_JONES_HPP_INCLUDED
#include <cmath>
#include <vec3.hpp>
#include <iostream>

/*Efficiency idea: replace ' r < cutoff_ ' with ' std::signbit(r - cutoff_) '. This avoids branching at the cost of many more instructions.*/


template <typename real_t>
class Lennard_Jones_Potential{
public:
  Lennard_Jones_Potential() : epsilon_(1), sigma_(1), cutoff_(2.5), derivative_offset_(derivative_(2.5)), offset_(clean_eval_(2.5)) {}
  Lennard_Jones_Potential(real_t epsilon, real_t sigma): epsilon_(epsilon), sigma_(sigma), cutoff_(2.5 * sigma), derivative_offset_(derivative_(2.5 * sigma)), offset_(clean_eval_(2.5 * sigma)) {}

  real_t operator()(real_t r){
    return (r < cutoff_) * (real_t(4) * epsilon_ * (std::pow(sigma_ / r, 12) - std::pow(sigma_ / r, 6)) - offset_ - derivative_offset_ * r);
  }
private:
  real_t derivative_(real_t r) {
    return 24 * epsilon_ * (-2 * std::pow(sigma_, 12) * std::pow(r, -13) + std::pow(sigma_, 6) * std::pow(r, -7));
  }
  real_t clean_eval_(real_t r) {
    return 4 * epsilon_ * (std::pow(sigma_, 12) * std::pow(r, -12) - std::pow(sigma_, 6) * std::pow(r, -6)) - derivative_offset_ * r;
  }
  const real_t epsilon_, sigma_;
  const real_t cutoff_;
  const real_t derivative_offset_, offset_;
};

template <typename real_t>
class Lennard_Jones_Force{
public:
  Lennard_Jones_Force(): epsilon_(1), sigma_(1), cutoff_(2.5), derivative_offset_(derivative_(2.5)) {}
  Lennard_Jones_Force(real_t epsilon, real_t sigma): epsilon_(epsilon), sigma_(sigma), cutoff_(2.5 * sigma), derivative_offset_(derivative_(2.5 * sigma)) {}

  Vec3<real_t> operator()(const Vec3<real_t>& rvec){ //force on i, rvec = r_j - r_i
    real_t r = rvec.r();
    return real_t((r < cutoff_) * (24 * epsilon_ * (-2 * std::pow(sigma_ / r, 14) * (1. / (sigma_ * sigma_)) + std::pow(sigma_ / r, 8) * (1. / (sigma_ * sigma_))) - derivative_offset_ / r)) * rvec;
  }
private:
  real_t derivative_(real_t r) {
    return 24 * epsilon_ * (-2 * std::pow(sigma_, 12) * std::pow(r, -13) + std::pow(sigma_, 6) * std::pow(r, -7));
  }
  const real_t epsilon_, sigma_;
  const real_t cutoff_;
  const real_t derivative_offset_;
};

#endif
