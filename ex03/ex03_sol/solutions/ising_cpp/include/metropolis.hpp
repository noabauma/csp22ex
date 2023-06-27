#pragma once

#include <array>

#include "ising_lattice.hpp"
#include "typedefs.hpp"

class Metropolis final : public IsingLattice {
public:
  Metropolis(int L, int dimension);

  void doSweep(int sweep_size) override;
  void setBeta(Real beta) override;
  void energize(int energy) override; 

private:
  void doStep();
  Real computeProb(int halo) const;
  void increaseEnergy(); 

  int sz = 2 * dimension_ + 1;
  std::vector<Real> exp_table_ = std::vector<Real>(sz);
};
