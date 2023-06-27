#pragma once

#include <array>

#include "ising_lattice.hpp"
#include "typedefs.hpp"

class Metropolis final : public IsingLattice {
public:
  Metropolis(int L);

  void doSweep(int sweep_size) override;

  void setBeta(Real beta) override;

private:
  void doStep();
  Real computeProb(int halo) const;

  // Note: lookup table size is 5 in 2-dimensions.
  std::array<Real, 5> exp_table_; 
};
