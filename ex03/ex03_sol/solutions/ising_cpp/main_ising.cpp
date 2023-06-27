#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>

#include "json11.hpp"
#include "metropolis.hpp"
#include "creutz.hpp"

int main(int argc, char **argv)
{
  const std::string inputname = argc > 1 ? argv[1] : "input.json";
  std::cout << "Reading " << inputname << "\n";
  json11::Json reader = json11::Json::parseFile(inputname);

  const std::string algorithm = reader["algorithm"].string_value();

  const int dimension = reader["dimension"].int_value(); // dimensionality of the lattice

  const int sweep_factor = reader["sweep-factor"].int_value(); // sweep_size is sweep_factor*L**dim.
  const int thermalization_factor = reader["thermalization-factor"].int_value();
  const int measurements = reader["measurements"].int_value(); // number of measurements.

  const int E_d_init = reader["initial E_d"].int_value(); // initial demon energy
  const int E_max = reader["E_max"].int_value(); // energy capacity of demon in Creutz algorithm.
  const int num_E = reader["num_energies"].int_value(); // energy capacity of demon in Creutz algorithm.

  // print current simulation method and dimension.
  std::cout << "Simulating " << dimension << "-d Ising model using the algorithm: " << algorithm << std::endl; 

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
    std::cout << std::endl;
    std::cout << "Size: " << L << std::endl; // print current linear system size

    if (algorithm == "Metropolis")
    {
      std::ofstream out("outputs/thdyn_" + std::to_string(dimension) + "d_L" + std::to_string(L) + ".txt");
      out << "# beta\tE\tM\tE2\tM2\n";

      std::ofstream e_file("outputs/energies_" + std::to_string(dimension) + "d_L" + std::to_string(L) + ".txt");
      e_file << "# beta\tE\n";
      e_file.precision(10);

      std::ofstream m_file("outputs/magnetizations_" + std::to_string(dimension) + "d_L" + std::to_string(L) +
                           ".txt");
      m_file << "# beta\tM\n";
      m_file.precision(10);

      Metropolis lattice(L, dimension);

      int id = 0;
      for (Real beta : betas)
      {
        std::cout << "Inverse temperature: " << beta << std::endl; // print current simulation temperature.
        lattice.setBeta(beta);

        e_file << beta << "\t";
        m_file << beta << "\t";

        // Thermalize.
        for (int i = 0; i < thermalization_factor; ++i)
        {
          lattice.doSweep(sweep_factor * lattice.n_);
        }
        lattice.markThermalized();

        // Measure.
        Real E(0), E2(0), M(0), M2(0);
        for (int i = 0; i < measurements; ++i)
        {
          lattice.doSweep(sweep_factor * lattice.n_);
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
        M /= measurements; // Store magnetisation density.
        M2 /= measurements;
        out << beta << "\t" << E << "\t" << M << "\t" << E2 << "\t" << M2
            << std::endl;
        ++id;
      }
    }

    else if (algorithm == "Creutz")
    {
      Creutz lattice(L, dimension);
      lattice.setEmax(E_max);

      int del_E = dimension * lattice.n_ / num_E;
      std::vector<int> energies(num_E-1);
      std::generate(energies.begin(), energies.end(), [e = -1, &del_E]() mutable { return del_E * e--; });

      std::ofstream out("outputs/thdyn_" + std::to_string(dimension) + "d_L" + std::to_string(L) + ".txt");
      out << "# energy\tEd\tE\tM\n";

      std::ofstream Ed_file("outputs/demon_energies_" + std::to_string(dimension) + "d_L" + std::to_string(L) + ".txt");
      Ed_file << "# energy\tEd\n";
      Ed_file.precision(10);

      std::ofstream e_file("outputs/energies_" + std::to_string(dimension) + "d_L" + std::to_string(L) + ".txt");
      e_file << "# energy\tE\n";
      e_file.precision(10);

      std::ofstream m_file("outputs/magnetizations_" + std::to_string(dimension) + "d_L" + std::to_string(L) + ".txt");
      m_file << "# energy\tM\n";
      m_file.precision(10);

      

      int id = 0;
      for (int energy : energies)
      {
        std::cout << "Energy: " << energy << std::endl; // print current simulation energy.

        Ed_file << energy << "\t";
        e_file << energy << "\t";
        m_file << energy << "\t";
        
        lattice.energize(energy); // Equilibrate by incerasing the energy.
        lattice.setEd(E_d_init); // demon with energy E_d_init is released to the system.

        // Measure.
        Real Ed(0); 
        Real E(0), M(0);
        for (int i = 0; i < measurements; ++i)
        {
          lattice.doSweep(sweep_factor * lattice.n_);

          Ed_file << lattice.getEd() << "\t";
          e_file << lattice.getE() << "\t";
          m_file << lattice.getM() << "\t";

          Ed += Real(lattice.getEd());
          E += lattice.getE();
          M += abs(lattice.getM());
        }
        Ed_file << "\n";
        e_file << "\n";
        m_file << "\n";

        Ed /= measurements;
        E /= measurements; // Store energy density.
        M /= measurements; // Store magnetisation density.

        // print current simulation temperature:
        Real T_eff = 4.0 / log(1.0 + 4.0 / Ed);
        std::cout << "Mean demon energy is: " << Ed << std::endl;
        std::cout << "Effective temperature is: " << T_eff << std::endl;

        out << energy << "\t" << Ed << "\t" << E << "\t" << M << std::endl;
        ++id;
      }
    }
  }
  return 0;
}
