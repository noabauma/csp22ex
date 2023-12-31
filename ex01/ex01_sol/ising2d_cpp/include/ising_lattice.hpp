// Base lattice class.
#pragma once

#include <array>
#include <random>
#include <vector>

#include "typedefs.hpp"

class IsingLattice
{
public:
  IsingLattice(int L);

  virtual void doSweep(int sweep_size) = 0;
  virtual void setBeta(Real beta) = 0;
  virtual void markThermalized() {}

  Real getE() const { return Real(E_) / n_; } // Return energy density.
  Real getM() const { return Real(M_) / n_; } // Return magnetization density.
  std::vector<std::int8_t> getX() const { return std::vector<std::int8_t>(spins_); }
  auto size() const { return n_; }

protected:
  int haloMagnetization(int idx) const;
  // Compute the halo magnetization only for sites with higher linear index.
  int rightMagnetization(int i, int j) const;
  void computeEandM();

  const std::string lattice_;
  const int L_;
  const int n_;
  Real beta_;
  std::vector<std::int8_t> spins_;

  std::mt19937_64 rng_;

  long int E_; // in units of J.
  long int M_; // not normalized.
};
