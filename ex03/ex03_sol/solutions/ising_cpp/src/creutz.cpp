#include "creutz.hpp"
#include <cassert>

Creutz::Creutz(int L, int dimension) : IsingLattice(L, dimension) {}

void Creutz::doSweep(const int sweep_size)
{
    for (int i = 0; i < sweep_size; ++i)
        doStep();
}

void Creutz::setBeta(Real beta) { beta_ = beta; } // TODO: remove this

void Creutz::increaseEnergy()
{
    std::uniform_int_distribution<int> distro_int(0, n_ - 1);
    const int candidate = distro_int(rng_);

    const int halo = haloMagnetization(candidate);
    const int s_old = spins_[candidate];
    const int delta_E = 2 * s_old * halo;

    if (delta_E > 0)
    {
        spins_[candidate] *= -1; // flip spin.
        M_ += -2 * s_old;
        E_ += delta_E;
    }
}

void Creutz::energize(int energy)
{
    for (int i = 0; i < n_; ++i)
    {
        spins_[i] = 1; // start with uniform initial configuration.
    }
    E_ = -dimension_ * n_;
    M_ = n_;
    while (energy > E_)
    {
        increaseEnergy();
    }
}

void Creutz::doStep()
{
    std::uniform_int_distribution<int> distro_int(0, n_ - 1);
    const int candidate = distro_int(rng_);

    const int halo = haloMagnetization(candidate);
    const int s_old = spins_[candidate];
    const int delta_E = 2 * s_old * halo;

    if (E_max_ >= Ed_ - delta_E and Ed_ >= delta_E)
    {
        spins_[candidate] *= -1; // flip spin.
        M_ += -2 * s_old;
        E_ += delta_E;
        Ed_ -= delta_E;
    }
}
