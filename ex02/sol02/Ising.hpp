/* Ising.hpp 2021
 * Author: Pascal Engeler
 * engelerp@phys.ethz.ch
 */

#ifndef ISING_HPP_INCLUDED
#define ISING_HPP_INCLUDED
#include <vector> //std::vector
#include <string> //std::string
#include <random> //std::uniform_real_distribution

//stepper types
using MRT_type = double;


/** Class to represent an Ising system.
  * Template arguments:
  * unsigned ISING_DIM -- Dimension of the system. Valid values: 2, 3.
  *
  * Workflow:
  * The class is designed to be used with the following workflow:
  *   1. Construct the system with the desired size, coupling and dimension
  *   2. Seed the RNG
  *   3. If measurements for different temperatures are desired, then for each temperature t:
  *     - Set up the measurement parameters
  *     - run the measurement / correlation analysis with the desired stepper
  *   4. Once all runs are finished, save the data / correlation data
**/
template<unsigned ISING_DIM>
class Ising{
public:

  /* Constructor
   * N -- Linear system size. Only powers of two are allowed (see Ising::step_(unsigned)).
   * J -- Coupling constant
   */
  Ising(unsigned N=8, double J=1.);
  Ising()=delete;

  /* Seed the random number generator
   * s -- Desired seed
   */
  void seed(int s);

  /* Set the parameters to be used for the next measurement
   * temperature -- System temperature
   * thermalize_steps -- Number of updates to be performed before measurements start
   * num_samples -- Number of samples to be collected
   * sample_stride -- Number of configurations discarded between successive samples
   */
  void set_measurement_parameters(double temperature, unsigned thermalize_steps, unsigned num_samples, unsigned sample_stride);

  /* Run measurement with previously specified parameters
   * Template arguments:
   * T -- specifies the stepper type to be used. Valid values: MRT_type for M(RT)^2
   */
  template <typename T>
  void run_measurement();

  /* Run correlation analysis with previously specified parameters.
   * This method measures after each update, sample_stride is ignored.
   * Template arguments:
   * T -- specifies the stepper type to be used. Valid values: MRT_type for M(RT)^2
   */
  template <typename T>
  void run_correlation_analysis();

  /* Save the data taken in the previous measurement as comma separated values.
   * The data is best analyzed using pandas.
   * filename -- File to which the data is saved.
   */
  void save_data(std::string filename) const;

  /* Save the data taken in the previous correlation analysis as comma separated values.
   * The data is best analyzed using pandas.
   * filename -- File to which the data is saved.
   */
  void save_correlation_data(std::string filename) const;


private:
  //private member functions
  void clear_samples_();
  void reserve_samples_();
  template<typename T>
  void step_(int num=1);
  void fetch_sample_();
  void extract_data_();
  //observable extraction
  double compute_energy_() const;
  double compute_magnetization_() const;
  double compute_magnetization_correlation_() const;
  //tracked observable initialization
  void initialize_energy_tracking_();
  void initialize_magnetization_tracking_();
  //reset spin configuration
  void reset_configuration_();

  //System Parameters
  const unsigned N_; //linear system size
  const double J_; //coupling
  std::vector<short int> spins_; //sites
  std::vector<double> exp_table_;
  //std::mt19937_64 rng_; //random number generator
  std::minstd_rand0 rng_;
  std::uniform_real_distribution<float> disd_; //double distribution
  std::uniform_int_distribution<int> disi_;

  //Sampling Parameters
  double temperature_; //system temperature
  unsigned thermalize_steps_, num_samples_, sample_stride_;
  long int tracked_magnetization_;
  double tracked_energy_;
  std::vector<double> energies_; //energies
  std::vector<double> sq_energies_; //squared energies
  std::vector<double> magnetizations_; //magnetizations
  std::vector<double> sq_magnetizations_; //squared magnetizations

  //Results
  std::vector<double> res_temperatures_;
  std::vector<double> res_energies_;
  std::vector<double> res_magnetizations_;
  std::vector<double> res_susceptibilities_;
  std::vector<double> res_capacities_;
  //Correlation Data
  std::vector<double> corr_temperatures_;
  std::vector<std::vector<double> > corr_nonlinear_magnetizations_;
  std::vector<std::vector<double> > corr_linear_phis_;
  std::vector<short int> corr_t0_configuration_;
  double corr_t0_magnetization_;
};

#endif
