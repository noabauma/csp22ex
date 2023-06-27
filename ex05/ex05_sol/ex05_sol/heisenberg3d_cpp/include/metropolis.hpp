#pragma once

#include <array>

#include "heisenberg_lattice.hpp"
#include "typedefs.hpp"

class Metropolis final : public HeisenbergLattice {
public:
  Metropolis(int L);
  void doSweep() override;

private:
  void doStep();

  std::uniform_int_distribution<int> distro_int_;
  std::uniform_real_distribution<Real> distro_real_;
};
