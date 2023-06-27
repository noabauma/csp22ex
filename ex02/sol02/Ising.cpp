/* Ising.cpp 2021
 * Author: Pascal Engeler
 * engelerp@phys.ethz.ch
 */

#include "Ising.hpp"
#include <fstream> //std::fstream
#include <string> //std::string
#include <vector>
#include <iostream> //std::cout
#include <cmath> //std::exp
#include <algorithm> //std::for_each
#include <bit> //std::popcount

template<unsigned ISING_DIM>
Ising<ISING_DIM>::Ising(unsigned N, double J): N_(N),
                                        J_(J),
                                        rng_(),
                                        disd_(0.f,1.f),
                                        disi_(0,N-1),
                                        temperature_(0.),
                                        thermalize_steps_(0),
                                        num_samples_(0),
                                        sample_stride_(0)
{
  //check if system length is a power of two, i.e. 1 bit is set
  if(std::popcount(N) != 1){
    throw "SizeException: System size is not a power of 2.";
  }
  if constexpr(ISING_DIM == 2u){
    spins_ = std::vector<short int>(N*N, -1);
    exp_table_ = std::vector<double>(9);
  }
  else if constexpr(ISING_DIM == 3u){
    spins_ = std::vector<short int>(N*N*N, -1);
    exp_table_ = std::vector<double>(13);
  }
  else{
    throw "DimensionException: System dimension invalid.";
  }

  initialize_energy_tracking_();
  initialize_magnetization_tracking_();
}

template<unsigned ISING_DIM>
void Ising<ISING_DIM>::seed(int s){
  rng_.seed(s);
}

template<unsigned ISING_DIM>
void Ising<ISING_DIM>::set_measurement_parameters(double temperature, unsigned thermalize_steps, unsigned num_samples, unsigned sample_stride){
  temperature_ = temperature;
  thermalize_steps_ = thermalize_steps;
  num_samples_ = num_samples;
  sample_stride_ = sample_stride;
  //fill in the exp table: exp_table_[s*h+4] = min(1., std::exp(-2*J*s*h/temperature))
  if constexpr(ISING_DIM == 2u){
    exp_table_[0] = 1.;
    exp_table_[1] = 1.;
    exp_table_[2] = 1.;
    exp_table_[3] = 1.;
    exp_table_[4] = 1.;
    exp_table_[6] = std::exp(-2.*J_*2./temperature);
    exp_table_[8] = std::exp(-2.*J_*4./temperature);
  }
  else if constexpr(ISING_DIM == 3u){
    exp_table_[0] = 1.;
    exp_table_[1] = 1.;
    exp_table_[2] = 1.;
    exp_table_[3] = 1.;
    exp_table_[4] = 1.;
    exp_table_[5] = 1.;
    exp_table_[6] = 1.;
    exp_table_[8] = std::exp(-2.*J_*2./temperature);
    exp_table_[10] = std::exp(-2.*J_*4./temperature);
    exp_table_[12] = std::exp(-2.*J_*6./temperature);
  }
}

template<unsigned ISING_DIM>
template<typename T>
void Ising<ISING_DIM>::run_measurement(){
  //prepare for new run
  clear_samples_();
  reserve_samples_();
  //thermalize
  for(int i = 0;  i < thermalize_steps_; ++i){
    step_<T>();
  }
  //measure
  for(int i = 0; i < num_samples_; ++i){
    //decorrelate
    for(int j = 0; j < sample_stride_; ++j){
      step_<T>();
    }
    //fetch measurement result
    fetch_sample_();
  }
  //get measurement result
  extract_data_();
}

template<unsigned ISING_DIM>
template<typename T>
void Ising<ISING_DIM>::run_correlation_analysis(){
  //reserve space for data
  corr_nonlinear_magnetizations_.push_back(std::vector<double>(thermalize_steps_, 0.));
  //corr_nonlinear_magnetizations_.back().reserve(thermalize_steps_);
  //nonlinear relaxation analysis
  //average over 50 markov chains
  int num_chains = 50;
  for(int k = 0; k < num_chains; ++k){
    //start from a fresh configuration
    reset_configuration_();
    for(int i = 0; i < thermalize_steps_; ++i){
      step_<T>();
      corr_nonlinear_magnetizations_.back()[i] += compute_magnetization_();
    }
  }
  auto start = corr_nonlinear_magnetizations_.back().begin();
  auto end = corr_nonlinear_magnetizations_.back().end();
  std::transform(start, end, start, [&](auto m){return m / num_chains;});
  //linear relaxation analysis
  //reserve space for data
  corr_linear_phis_.push_back(std::vector<double>());
  corr_linear_phis_.back().reserve(num_samples_);
  //store the reference configuration ('t0')
  corr_t0_configuration_ = spins_;
  corr_t0_magnetization_ = compute_magnetization_();
  double denominator = 1. - corr_t0_magnetization_*corr_t0_magnetization_;

  for(int i = 0; i < num_samples_; ++i){
    step_<T>();
    double phi = (compute_magnetization_correlation_() - corr_t0_magnetization_*corr_t0_magnetization_) / denominator;
    corr_linear_phis_.back().push_back(phi);
  }
  corr_temperatures_.push_back(temperature_);
}

template<unsigned ISING_DIM>
void Ising<ISING_DIM>::save_data(std::string filename) const{
  std::ofstream file (filename);
  if(!file.is_open()){
    std::cout << "FileError: Failed to open file " << filename << "!" << std::endl;
  }
  else{
    //Header
    file << "Temperature,Energy,Magnetization,Susceptibility,Capacity,Dimension,Length" << std::endl;
    //data
    //on the first line, we save parameters
    file  << res_temperatures_[0] << ","
          << res_energies_[0] << ","
          << res_magnetizations_[0] << ","
          << res_susceptibilities_[0] << ","
          << res_capacities_[0] << ","
          << ISING_DIM << ","
          << N_
          << std::endl;
    for(int i = 1; i < res_temperatures_.size(); ++i){
      file  << res_temperatures_[i] << ","
            << res_energies_[i] << ","
            << res_magnetizations_[i] << ","
            << res_susceptibilities_[i] << ","
            << res_capacities_[i]
            << std::endl;
    }
    file.close();
  }
}

template<unsigned ISING_DIM>
void Ising<ISING_DIM>::save_correlation_data(std::string filename) const{
  std::ofstream file(filename);
  if(!file.is_open()){
    std::cout << "FileError: Failed to open file " << filename << "!" << std::endl;
  }
  else{
    //Header
    file << "Number Of Temperatures,";
    for(int i = 0; i < corr_temperatures_.size(); ++i){
      file << "Temperature " << i << ",Thermalization Magnetization " << i << ",Linear Phi " << i;
      if(i != corr_temperatures_.size()-1)
        file << ",";
    }
    file << std::endl;
    //First data line contains temperatures
    file << corr_temperatures_.size() << ",";
    for(int i = 0; i < corr_temperatures_.size(); ++i){
      file  << corr_temperatures_[i] << ","
            << corr_nonlinear_magnetizations_[i][0] << ","
            << corr_linear_phis_[i][0];
            if(i != corr_temperatures_.size()-1)
              file << ",";
    }
    file << std::endl;
    //Data
    for(int i = 0; i < std::max(corr_nonlinear_magnetizations_[0].size(), corr_linear_phis_[0].size()); ++i){
      file << ',';
      for(int j = 0; j < corr_temperatures_.size(); ++j){
        file << ",";
        if(i < corr_nonlinear_magnetizations_[j].size())
          file << corr_nonlinear_magnetizations_[j][i];
        file << ",";
        if(i < corr_linear_phis_[j].size())
          file << corr_linear_phis_[j][i];
        if(j < corr_temperatures_.size()-1)
          file << ",";
      }
      file << std::endl;
    }
    file.close();
  }
}

template<unsigned ISING_DIM>
void Ising<ISING_DIM>::clear_samples_(){
  energies_.clear();
  sq_energies_.clear();
  magnetizations_.clear();
  sq_magnetizations_.clear();
}

template<unsigned ISING_DIM>
void Ising<ISING_DIM>::reserve_samples_(){
  energies_.reserve(num_samples_);
  sq_energies_.reserve(num_samples_);
  magnetizations_.reserve(num_samples_);
  sq_magnetizations_.reserve(num_samples_);
}

template<unsigned ISING_DIM>
template<typename T>
void Ising<ISING_DIM>::step_(int num){
  //This is where the dimension limitation to powers of 2 comes from.
  //The restriction is used to implement an efficient modulo.
  //TODO: Lift this restriction by making the dimension a template argument, and specializing the neighbour calculation on it.
  if constexpr(std::is_same_v<T, MRT_type> && ISING_DIM == 2){
    //find spin to update
    int ix = disi_(rng_);
    int iy = disi_(rng_);
    //calculate neighboursum with efficient modulo
    int h = spins_[iy*N_+((ix-1+N_)&(N_-1))] + spins_[iy*N_+((ix+1)&(N_-1))] +  //along x
            spins_[((iy-1+N_)&(N_-1))*N_+ix] + spins_[((iy+1)&(N_-1))*N_+ix];   //along y
    //update spin
    short int old_spin = spins_[iy*N_+ix];
    spins_[iy*N_+ix] *= static_cast<short int>( exp_table_[h*spins_[iy*N_+ix]+4] + disd_(rng_) )*(-2) + 1;
    short int spin_delta = spins_[iy*N_+ix] - old_spin; //this is either 0, or +- 2
    tracked_magnetization_ += spin_delta;
    tracked_energy_ += static_cast<double>(std::abs(spin_delta) * J_ * h * old_spin);
  }
  else if constexpr(std::is_same_v<T, MRT_type> && ISING_DIM == 3){
    int ix = disi_(rng_);
    int iy = disi_(rng_);
    int iz = disi_(rng_);
    //selected spin is spins_[iz*N_*N_ + iy*N_ + ix]
    //neighbour field
    int h = spins_[iz*N_*N_ + iy*N_ + ((ix+1)&(N_-1))] + spins_[iz*N_*N_ + iy*N_ + ((ix-1+N_)&(N_-1))] +       //along x
            spins_[iz*N_*N_ + ((iy+1)&(N_-1))*N_ + ix] + spins_[iz*N_*N_ + ((iy-1+N_)&(N_-1))*N_ + ix] + //along y
            spins_[((iz+1)&(N_-1))*N_*N_ + iy*N_ + ix] + spins_[((iz-1+N_)&(N_-1))*N_*N_ + iy*N_ + ix];  //along z
    //update spin
    short int old_spin = spins_[iz*N_*N_ + iy*N_ + ix];
    spins_[iz*N_*N_ + iy*N_ + ix] *= static_cast<short int>( exp_table_[h*spins_[iz*N_*N_ + iy*N_ + ix]+6] + disd_(rng_) )*(-2) + 1;
    short int spin_delta = spins_[iz*N_*N_ + iy*N_ + ix] - old_spin; //this is either 0, or +- 2
    tracked_magnetization_ += spin_delta;
    tracked_energy_ += static_cast<double>(std::abs(spin_delta) * J_ * h * old_spin);
  }
  else{
    throw "Ising::step_: Invalid Stepper-Dimension Combination Selected.";
  }
}

template<unsigned ISING_DIM>
double Ising<ISING_DIM>::compute_energy_() const{
  /*
  double energy = 0;
  for(int ix = 0; ix < N_; ++ix){
    for(int iy = 0; iy < N_; ++iy){
      if constexpr(ISING_DIM == 2){
        energy -= J_*spins_[iy*N_+ix]*( spins_[iy*N_+((ix-1+N_)&(N_-1))] + spins_[iy*N_+((ix+1)&(N_-1))] +
                                        spins_[((iy-1+N_)&(N_-1))*N_+ix] + spins_[((iy+1)&(N_-1))*N_+ix]);
      }
      else if constexpr(ISING_DIM == 3){
        for(int iz = 0; iz < N_; ++iz){
          energy -= J_*spins_[iz*N_*N_+iy*N_+ix]*(spins_[iz*N_*N_ + iy*N_ + ((ix+1)&(N_-1))] + spins_[iz*N_*N_ + iy*N_ + ((ix-1+N_)&(N_-1))] +
                                                  spins_[iz*N_*N_ + ((iy+1)&(N_-1))*N_ + ix] + spins_[iz*N_*N_ + ((iy-1+N_)&(N_-1))*N_ + ix] +
                                                  spins_[((iz+1)&(N_-1))*N_*N_ + iy*N_ + ix] + spins_[((iz-1+N_)&(N_-1))*N_*N_ + iy*N_ + ix]);
        }
      }
    }
  }
  return energy/(2.*static_cast<double>(spins_.size()));
  */
  //TODO: This works as long as we can track the energy. Else, specialize it.
  return tracked_energy_ / static_cast<double>(spins_.size());
}

template<unsigned ISING_DIM>
double Ising<ISING_DIM>::compute_magnetization_() const{
  /*
  double magnetization = 0;
  std::for_each(spins_.begin(), spins_.end(), [&magnetization](short int s){magnetization += static_cast<double>(s);});
  return magnetization/static_cast<double>(spins_.size());
  */
  //TODO: This works as long as we can track the magnetization. Else, specialize it.
  return static_cast<double>(tracked_magnetization_) / static_cast<double>(spins_.size());
}

template<unsigned ISING_DIM>
void Ising<ISING_DIM>::reset_configuration_(){
  if constexpr(ISING_DIM == 2u){
    spins_ = std::vector<short int>(N_*N_, -1);
  }
  else if constexpr(ISING_DIM == 3u){
    spins_ = std::vector<short int>(N_*N_*N_, -1);
  }
  initialize_energy_tracking_();
  initialize_magnetization_tracking_();
}

template<unsigned ISING_DIM>
double Ising<ISING_DIM>::compute_magnetization_correlation_() const{
  double result = 0.;
  int j = 0;
  std::for_each(spins_.begin(), spins_.end(), [&result, &j, this](short int s){result += static_cast<double>(s*corr_t0_configuration_[j++]);});
  return result/static_cast<double>(spins_.size());
}

template<unsigned ISING_DIM>
void Ising<ISING_DIM>::initialize_energy_tracking_(){
  tracked_energy_ = 0.;
  for(int ix = 0; ix < N_; ++ix){
    for(int iy = 0; iy < N_; ++iy){
      if constexpr(ISING_DIM == 2){
        tracked_energy_ -= J_*spins_[iy*N_+ix]*( spins_[iy*N_+((ix-1+N_)&(N_-1))] + spins_[iy*N_+((ix+1)&(N_-1))] +
                                        spins_[((iy-1+N_)&(N_-1))*N_+ix] + spins_[((iy+1)&(N_-1))*N_+ix]);
      }
      else if constexpr(ISING_DIM == 3){
        for(int iz = 0; iz < N_; ++iz){
          tracked_energy_ -= J_*spins_[iz*N_*N_+iy*N_+ix]*(spins_[iz*N_*N_ + iy*N_ + ((ix+1)&(N_-1))] + spins_[iz*N_*N_ + iy*N_ + ((ix-1+N_)&(N_-1))] +
                                                  spins_[iz*N_*N_ + ((iy+1)&(N_-1))*N_ + ix] + spins_[iz*N_*N_ + ((iy-1+N_)&(N_-1))*N_ + ix] +
                                                  spins_[((iz+1)&(N_-1))*N_*N_ + iy*N_ + ix] + spins_[((iz-1+N_)&(N_-1))*N_*N_ + iy*N_ + ix]);
        }
      }
    }
  }
  tracked_energy_ /= 2.;
}

template<unsigned ISING_DIM>
void Ising<ISING_DIM>::initialize_magnetization_tracking_(){
  double magnetization = 0.;
  std::for_each(spins_.begin(), spins_.end(), [&magnetization](short int s){magnetization += static_cast<double>(s);});
  tracked_magnetization_ = magnetization;
}

template<unsigned ISING_DIM>
void Ising<ISING_DIM>::fetch_sample_(){
  double energy = compute_energy_();
  double magnetization = std::abs(compute_magnetization_());
  energies_.push_back(energy);
  sq_energies_.push_back(energy*energy);
  magnetizations_.push_back(magnetization);
  sq_magnetizations_.push_back(magnetization*magnetization);
}

template<unsigned ISING_DIM>
void Ising<ISING_DIM>::extract_data_(){
  double energy=0., sq_energy=0.;
  double magnetization=0., sq_magnetization=0.;
  std::for_each(energies_.begin(), energies_.end(), [&energy](double e){energy+=e;});
  energy /= static_cast<double>(energies_.size());
  std::for_each(sq_energies_.begin(), sq_energies_.end(), [&sq_energy](double esq){sq_energy+=esq;});
  sq_energy /= static_cast<double>(sq_energies_.size());
  std::for_each(magnetizations_.begin(), magnetizations_.end(), [&magnetization](double m){magnetization+=m;});
  magnetization /= static_cast<double>(magnetizations_.size());
  std::for_each(sq_magnetizations_.begin(), sq_magnetizations_.end(), [&sq_magnetization](double msq){sq_magnetization+=msq;});
  sq_magnetization /= static_cast<double>(magnetizations_.size());

  res_temperatures_.push_back(temperature_);
  res_energies_.push_back(energy);
  res_magnetizations_.push_back(magnetization);
  res_susceptibilities_.push_back((sq_magnetization-magnetization*magnetization)*static_cast<double>(spins_.size())/temperature_);
  res_capacities_.push_back((sq_energy-energy*energy)*static_cast<double>(spins_.size())/(temperature_*temperature_));
}



//force compiler to produce the desired template blueprints
template void Ising<2u>::run_measurement<MRT_type>();
template void Ising<3u>::run_measurement<MRT_type>();
template void Ising<2u>::run_correlation_analysis<MRT_type>();
template void Ising<3u>::run_correlation_analysis<MRT_type>();
template class Ising<2u>;
template class Ising<3u>;
