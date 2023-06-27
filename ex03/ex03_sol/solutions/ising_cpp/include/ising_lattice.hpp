// Base lattice class.
#pragma once

#include <array>
#include <random>
#include <vector>

#include "typedefs.hpp"

class IsingLattice
{
public:
  IsingLattice(int L, int dimension);

  virtual void doSweep(int sweep_size) = 0;
  virtual void setBeta(Real beta) = 0; // set the temperature of the canonical ensemble.
  virtual void energize(int energy) = 0;
  virtual void markThermalized() {}

  void setEd(int Ed){ Ed_ = Ed; } // set initial demon energy.
  void setEmax(int E_max) { E_max_ = E_max; } // set initial demon energy.

  Real getEd() const { return int(Ed_); } // Return demon energy for the Creutz algorithm.
  Real getE() const { return Real(E_) / n_; } // Return energy density.
  Real getM() const { return Real(M_) / n_; } // Return magnetization density.
  // get the full spin configuration as a flat list:
  std::vector<std::int8_t> getX() const { return std::vector<std::int8_t>(spins_); } 
  auto size() const { return n_; }

  const int n_; // system size

protected:
  int haloMagnetization(int idx) const;
  // Compute the halo magnetization only for sites with higher linear index:
  int rightMagnetization(int idx) const;
  void computeEandM();

  const std::string lattice_;
  const int dimension_;
  const int L_; // linear system size
  Real beta_; // temperature of the canonical ensemble

  std::vector<std::int8_t> spins_; // the full spin configuration.

  std::mt19937_64 rng_; // random number generator

  long int E_max_; // energy capacity of the demon in the Creutz algorithm.
  long int Ed_; // demon energy in the Creutz algorithm.
  long int E_; // in units of J.
  long int M_; // not normalized.
};
