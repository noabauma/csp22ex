#include "metropolis.hpp"

#include <cassert>

Metropolis::Metropolis(int L)
    : HeisenbergLattice(L), distro_int_(0, n_ - 1), distro_real_(0, 1) {}

void Metropolis::doSweep() {
  for (int i = 0; i < n_; ++i)
    doStep();
}

void Metropolis::doStep() {
  const int candidate = distro_int_(rng_);
  const Spin s_new = Spin::random(rng_); // TODO propose the new spin as a small deformation of the original one for higher acceptance probability!

  const Spin spin_diff = s_new - spins_[candidate];

  const Real delta_E = -spin_diff * haloMagnetization(candidate);
  const Real prob = std::exp(-delta_E * beta_);

  const bool accept = distro_real_(rng_) < prob;

  if (accept) {
    spins_[candidate] = s_new; // update spin.
    M_ = M_ + spin_diff;
    E_ += delta_E;
  }
}
