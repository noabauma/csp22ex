// Base lattice class.
#pragma once

#include <array>
#include <random>
#include <vector>

#include "spin.hpp"
#include "typedefs.hpp"

class HeisenbergLattice {
public:
  HeisenbergLattice(int L);

  virtual void doSweep() = 0;
  virtual void markThermalized() {}

  Real getE() const { return E_; }

  // Return magnetization norm.
  Real getM2() const { return std::sqrt(M_.norm2()); }

  auto size() const { return n_; }

  void setBeta(Real beta) { beta_ = beta; }

protected:
  Spin haloMagnetization(int idx) const;
  // Compute the halo magnetization only for sites with higher linear index.
  Spin rightMagnetization(int i, int j, int k) const;
  void computeEandM();

  const int L_;
  const int n_;
  Real beta_;
  std::vector<Spin> spins_;

  std::mt19937_64 rng_;

  Real E_; // in units of J.
  Spin M_; // not normalized.
};
