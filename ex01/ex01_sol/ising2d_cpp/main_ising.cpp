#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

#include "json11.hpp"
#include "metropolis.hpp"

int main(int argc, char **argv)
{
  const std::string inputname = argc > 1 ? argv[1] : "input.json";
  std::cout << "Reading " << inputname << "\n";
  json11::Json reader = json11::Json::parseFile(inputname);

  const int sweep_factor = reader["sweep-factor"].int_value(); // sweep_size is sweep_factor*L**2.
  const int thermalization_factor = reader["thermalization-factor"].int_value();
  const int measurements = reader["measurements"].int_value();

  std::vector<int> Ls;
  for (const auto &l : reader["Ls"].array_items())
  {
    Ls.push_back(l.int_value());
  }

  std::vector<Real> betas;
  for (const auto &T : reader["Ts"].array_items())
  {
    betas.push_back(1. / T.number_value());
  }
  std::sort(betas.begin(), betas.end());

  for (int L : Ls)
  {
    Metropolis lattice(L);
    std::cout << "Size: " << L << std::endl;

    std::ofstream out("outputs/thdyn_L" + std::to_string(L) + ".txt");
    out << "# beta\tE\tM\tE2\tM2\n";

    std::ofstream e_file("outputs/energies_L" + std::to_string(L) + ".txt");
    e_file << "# beta\tE\n";
    e_file.precision(10);

    std::ofstream m_file("outputs/magnetizations_L" + std::to_string(L) +
                         ".txt");
    m_file << "# beta\tM\n";
    m_file.precision(10);

    int id  = 0;
    for (Real beta : betas)
    {
      std::cout << "Inverse temperature: " << beta << std::endl;
      lattice.setBeta(beta);

      e_file << beta << "\t";
      m_file << beta << "\t";

      // Thermalize.
      for (int i = 0; i < thermalization_factor; ++i)
      {
        lattice.doSweep(sweep_factor * L * L);
      }

      lattice.markThermalized();

      // Measure.
      Real E(0), E2(0), M(0), M2(0);
      for (int i = 0; i < measurements; ++i)
      {
        lattice.doSweep(sweep_factor * L * L);
        e_file << lattice.getE() << "\t";
        m_file << lattice.getM() << "\t";

        E += lattice.getE();
        E2 += lattice.getE() * lattice.getE();
        M += abs(lattice.getM());
        M2 += lattice.getM() * lattice.getM();
      }
      e_file << "\n";
      m_file << "\n";

      E /= measurements; // Store energy density.
      E2 /= measurements;
      M /= measurements;
      M2 /= measurements;
      out << beta << "\t" << E << "\t" << M << "\t" << E2 << "\t" << M2
          << std::endl;
      ++id;
    }
  }

  return 0;
}
