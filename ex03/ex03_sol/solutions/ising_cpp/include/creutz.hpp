#pragma once

#include <array>

#include "ising_lattice.hpp"
#include "typedefs.hpp"

class Creutz final : public IsingLattice
{
public:
	Creutz(int L, int dimension);
	void doSweep(int sweep_size) override;
	void setBeta(Real beta) override; // TODO: remove this
	void energize(int energy) override;
	// void thermometer();

private:
	void doStep();
	void increaseEnergy();
	Real computeProb(int halo) const;	
};
