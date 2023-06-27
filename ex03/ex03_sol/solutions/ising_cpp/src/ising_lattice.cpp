#include "ising_lattice.hpp"
#include <cassert>
#include <cstdarg>

int intPow(int x, int p) // integer power operation
{
  int i = 1;
  for (int j = 1; j <= p; j++)
    i *= x;
  return i;
}

IsingLattice::IsingLattice(int L, int dimension)
    : L_(L), dimension_(dimension), n_(intPow(L, dimension)),
      spins_(n_), rng_(std::random_device()()), E_(0), M_(0)
{
  std::uniform_int_distribution<std::int8_t> distro(0, 1);
  for (int i = 0; i < n_; ++i)
  {
    spins_[i] = 2 * distro(rng_) - 1; // draw binary random variables
  }
  computeEandM();
}

void IsingLattice::computeEandM()
{
  E_ = M_ = 0;

  for (int linindex = 0; linindex < n_; ++linindex)
  {
    M_ += spins_[linindex];
    E_ -= spins_[linindex] * rightMagnetization(linindex);
  }
}

int IsingLattice::haloMagnetization(int idx) const
{
  // Measures the local magnetic field on a spin due to its nearest-neighbours.
  // Take care of periodic boundary conditions
  // Note: don't exclude that implementing this in terms of '%' operator is
  // faster on certain architectures.
  auto pbs = [&](const int i) { //TODO: use s[i - 1 + L & (L - 1)]
    if (i >= 0 && i < L_)
      return i;
    else if (i >= L_)
      return i - L_;
    else
      return i + L_;
  };

  auto index = [&](int num, ...) {
    va_list args;
    va_start(args, num);

    int id = 0;
    for (int d = 0; d < num; d++)
    {
      int i = va_arg(args, int);
      id += intPow(L_, d + 1) / L_ * pbs(i);
    }
    return id;
    va_end(args);
  };

  if (dimension_ == 2)
  {
    // Compute linear indices for 2-d lattice.
    const int j = idx / L_;
    const int i = idx - L_ * j;

    return spins_[index(dimension_, i + 1, j)] + spins_[index(dimension_, i - 1, j)] +
           spins_[index(dimension_, i, j + 1)] + spins_[index(dimension_, i, j - 1)];
  }
  else if (dimension_ == 3)
  {
    // Compute linear indices for 3-d lattice.
    const int k = idx / (L_ * L_);
    idx -= k * (L_ * L_);
    const int j = idx / L_;
    const int i = idx - L_ * j;

    return spins_[index(dimension_, i + 1, j, k)] + spins_[index(dimension_, i - 1, j, k)] +
           spins_[index(dimension_, i, j + 1, k)] + spins_[index(dimension_, i, j - 1, k)] +
           spins_[index(dimension_, i, j, k + 1)] + spins_[index(dimension_, i, j, k - 1)];
  }
  else { return spins_[idx + 1] + spins_[idx - 1]; } // Halo magnetisation for the 1-d chain
}

int IsingLattice::rightMagnetization(int idx) const
{
  auto pbs = [&](const int i) {
    if (i >= 0 && i < L_)
      return i;
    else if (i >= L_)
      return i - L_;
    else
      return i + L_;
  };

  auto index = [&](int num, ...) {
    va_list args;
    va_start(args, num);

    int id = 0;
    for (int d = 0; d < num; d++)
    {
      int i = va_arg(args, int);
      id += intPow(L_, d + 1) / L_ * pbs(i);
    }
    return id;
    va_end(args);
  };

  if (dimension_ == 2)
  {
    // Compute unrolled indices for the 2-d lattice.
    const int j = idx / L_;
    const int i = idx - L_ * j;

    return spins_[index(dimension_, i + 1, j)] + spins_[index(dimension_, i, j + 1)];
  }
  else if (dimension_ == 3)
  {
    // Compute unrolled indices for the 3-d lattice.
    const int k = idx / (L_ * L_);
    idx -= k * (L_ * L_);
    const int j = idx / L_;
    const int i = idx - L_ * j;

    return spins_[index(dimension_, i + 1, j, k)] +
           spins_[index(dimension_, i, j + 1, k)] +
           spins_[index(dimension_, i, j, k + 1)];
  }
  else{ return spins_[idx + 1]; } // Right magnetisation for the 1-d chain
}