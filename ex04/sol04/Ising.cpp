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
#include <chrono> //clocks

template<unsigned ISING_DIM>
Ising<ISING_DIM>::Ising(unsigned N, double J): N_(N),
                                        J_(J),
                                        rng_(),
                                        disd_(0.f,1.f),
                                        disi_(0,N-1),
                                        temperature_(0.),
                                        thermalize_steps_(0),
                                        num_samples_(0),
                                        sample_stride_(0),
                                        sw_label_(std::pow(N, ISING_DIM)),
                                        sw_labels_(std::pow(N, ISING_DIM)),
                                        sw_flipper_(std::pow(N, ISING_DIM))
{
  //check if system length is a power of two, i.e. exactly 1 bit is set
  if(std::popcount(N) != 1){
    throw "SizeException: System size is not a power of 2.";
  }
  if constexpr(ISING_DIM == 2u){
    spins_ = std::vector<short int>(N*N, -1);
    exp_table_ = std::vector<double>(9);
    cl_wolff_queue_.reserve(N*N); //TODO: That's very memory inefficient
  }
  else if constexpr(ISING_DIM == 3u){
    spins_ = std::vector<short int>(N*N*N, -1);
    exp_table_ = std::vector<double>(13);
    cl_wolff_queue_.reserve(N*N*N); //TODO: That's very memory inefficient
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

  //cluster growth probability for current temperature
  cluster_growth_probability_ = 1. - std::exp(-2.*J_/temperature);

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
    fetch_sample_<T>();
  }
  //get measurement result
  extract_data_();
}

template<unsigned ISING_DIM>
template<typename T>
void Ising<ISING_DIM>::run_correlation_analysis(){
  //start from a fresh configuration
  reset_configuration_();
  //reserve space for data
  corr_nonlinear_magnetizations_.push_back(std::vector<double>());
  corr_nonlinear_magnetizations_.back().reserve(thermalize_steps_);
  //nonlinear relaxation analysis
  for(int i = 0; i < thermalize_steps_; ++i){
    step_<T>();
    corr_nonlinear_magnetizations_.back().push_back(compute_magnetization_());
  }
  //linear relaxation analysis
  //reserve space for data
  corr_linear_phis_.push_back(std::vector<double>());
  corr_linear_phis_.back().reserve(num_samples_);
  corr_linear_energies_.push_back(std::vector<double>());
  corr_linear_energies_.back().reserve(num_samples_);
  corr_linear_magnetizations_.push_back(std::vector<double>());
  corr_linear_magnetizations_.back().reserve(num_samples_);
  //store the reference configuration ('t0')
  corr_t0_configuration_ = spins_;
  corr_t0_magnetization_ = compute_magnetization_();
  double denominator = 1. - corr_t0_magnetization_*corr_t0_magnetization_;

  if constexpr(std::is_same_v<T, MRT_type>){
    for(int i = 0; i < num_samples_; ++i){
      //we perform one sweep between measurements
      unsigned intermediate_updates = 0;
      while(intermediate_updates < N_*N_*N_){
        intermediate_updates += step_<T>();
      }
      double phi = (compute_magnetization_correlation_() - corr_t0_magnetization_*corr_t0_magnetization_) / denominator;
      corr_linear_phis_.back().push_back(phi);
      corr_linear_energies_.back().push_back(compute_energy_<T>());
      corr_linear_magnetizations_.back().push_back(compute_magnetization_());
    }
    corr_temperatures_.push_back(temperature_);
  }
  else if constexpr(std::is_same_v<T, Wolff_type>){
    //first try to find average number of sites updated per cluster step
    unsigned num_sites_updated = 0;
    for(int i = 0; i < 500; ++i){
      num_sites_updated += step_<T>();
    }
    double steps_per_sweep = static_cast<double>(N_*N_*N_) / (static_cast<double>(num_sites_updated) / 500.);
    std::cout << "steps per sweep: " << steps_per_sweep << std::endl;
    int sites_per_step = static_cast<int>(static_cast<double>(num_sites_updated)/500.);
    for(int i = 0; i < num_samples_; ++i){
      //we perform one sweep between measurements, using our best estimates
      int site_count = 0;
      for(int j = 0; (j < static_cast<int>(steps_per_sweep)) && (N_*N_*N_ > site_count + sites_per_step); ++j){
        site_count += step_<T>();
      }
      double phi = (compute_magnetization_correlation_() - corr_t0_magnetization_*corr_t0_magnetization_) / denominator;
      corr_linear_phis_.back().push_back(phi);
      corr_linear_energies_.back().push_back(compute_energy_<T>());
      corr_linear_magnetizations_.back().push_back(compute_magnetization_());
    }
    corr_temperatures_.push_back(temperature_);
  }
  else if constexpr(std::is_same_v<T, SwendsenWang_type>){
    //first try to find average number of sites updated per cluster step
    unsigned num_sites_updated = 0;
    for(int i = 0; i < 500; ++i){
      num_sites_updated += step_<T>();
    }
    double steps_per_sweep = static_cast<double>(N_*N_*N_) / (static_cast<double>(num_sites_updated) / 500.);
    std::cout << "steps per sweep: " << steps_per_sweep << std::endl;
    int sites_per_step = static_cast<int>(static_cast<double>(num_sites_updated)/500.);
    for(int i = 0; i < num_samples_; ++i){
      //we perform one sweep between measurements, using our best estimates
      int site_count = 0;
      for(int j = 0; (j < static_cast<int>(steps_per_sweep)) && (N_*N_*N_ > site_count + sites_per_step); ++j){
        site_count += step_<T>();
      }
      double phi = (compute_magnetization_correlation_() - corr_t0_magnetization_*corr_t0_magnetization_) / denominator;
      corr_linear_phis_.back().push_back(phi);
      corr_linear_energies_.back().push_back(compute_energy_<T>());
      corr_linear_magnetizations_.back().push_back(compute_magnetization_());
    }
    corr_temperatures_.push_back(temperature_);
  }
}

template<unsigned ISING_DIM>
template<typename T>
void Ising<ISING_DIM>::run_benchmark(){
  std::string method_name = "";
  if constexpr(std::is_same_v<T, MRT_type>){
    method_name = "MRT";
  }
  else if constexpr(std::is_same_v<T, Wolff_type>){
    method_name = "Wolff";
  }
  else if constexpr(std::is_same_v<T, SwendsenWang_type>){
    method_name = "SW";
  }

  for(int i = 0; i < thermalize_steps_; ++i){
    step_<T>();
  }

  unsigned long long spins_updated = 0;
  auto start = std::chrono::high_resolution_clock::now();
  for(int i = 0; i < num_samples_; ++i){
    spins_updated += step_<T>();
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> duration_ms = end - start;
  bench_data_.push_back({method_name, ISING_DIM, N_, temperature_, duration_ms.count(), spins_updated, num_samples_});
}

template<unsigned ISING_DIM>
void Ising<ISING_DIM>::save_data(std::string filename) const{
  std::ofstream file (filename);
  if(!file.is_open()){
    std::cout << "FileError: Failed to open file " << filename << "!" << std::endl;
  }
  else{
    //Header
    file << "Temperature,Energy,Magnetization,Susceptibility,Capacity,BinderCumulant,Dimension,Length" << std::endl;
    //data
    //on the first line, we save parameters
    file  << res_temperatures_[0] << ","
          << res_energies_[0] << ","
          << res_magnetizations_[0] << ","
          << res_susceptibilities_[0] << ","
          << res_capacities_[0] << ","
          << res_binder_[0] << ","
          << ISING_DIM << ","
          << N_
          << std::endl;
    for(int i = 1; i < res_temperatures_.size(); ++i){
      file  << res_temperatures_[i] << ","
            << res_energies_[i] << ","
            << res_magnetizations_[i] << ","
            << res_susceptibilities_[i] << ","
            << res_capacities_[i] << ","
            << res_binder_[i]
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
      file << "Temperature " << i << ",Thermalization Magnetization " << i << ",Linear Phi " << i << ",Linear Energy " << i << ",Linear Magnetization " << i;
      if(i != corr_temperatures_.size()-1)
        file << ",";
    }
    file << std::endl;
    //First data line contains temperatures
    file << corr_temperatures_.size() << ",";
    for(int i = 0; i < corr_temperatures_.size(); ++i){
      file  << corr_temperatures_[i] << ","
            << corr_nonlinear_magnetizations_[i][0] << ","
            << corr_linear_phis_[i][0] << ","
            << corr_linear_energies_[i][0] << ","
            << corr_linear_magnetizations_[i][0];
            if(i != corr_temperatures_.size()-1)
              file << ",";
    }
    file << std::endl;
    //Data
    for(int i = 1; i < std::max(corr_nonlinear_magnetizations_[0].size(), corr_linear_phis_[0].size()); ++i){
      file << ',';
      for(int j = 0; j < corr_temperatures_.size(); ++j){
        file << ",";
        if(i < corr_nonlinear_magnetizations_[j].size())
          file << corr_nonlinear_magnetizations_[j][i];
        file << ",";
        if(i < corr_linear_phis_[j].size())
          file << corr_linear_phis_[j][i];
        file << ",";
        if(i < corr_linear_energies_[j].size())
          file << corr_linear_energies_[j][i];
        file << ",";
        if(i < corr_linear_magnetizations_[j].size())
          file << corr_linear_magnetizations_[j][i];
        if(j < corr_temperatures_.size()-1)
          file << ",";
      }
      file << std::endl;
    }
    file.close();
  }
}

template<unsigned ISING_DIM>
void Ising<ISING_DIM>::save_benchmark_data(std::string filename) const{
  std::ofstream file(filename);
  if(!file.is_open()){
    std::cout << "FileError: Failed to open file " << filename << "!" << std::endl;
  }
  else{
    //Header
    file << "Method,Dimension,Linear System Size,Temperature,Runtime,Number Sites Updated,Number Updates Performed\n";
    for(auto item: bench_data_){
      file  << std::get<0>(item) << ","
            << std::get<1>(item) << ","
            << std::get<2>(item) << ","
            << std::get<3>(item) << ","
            << std::get<4>(item) << ","
            << std::get<5>(item) << ","
            << std::get<6>(item) << "\n";
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
  qt_magnetizations_.clear();
}

template<unsigned ISING_DIM>
void Ising<ISING_DIM>::reserve_samples_(){
  energies_.reserve(num_samples_);
  sq_energies_.reserve(num_samples_);
  magnetizations_.reserve(num_samples_);
  sq_magnetizations_.reserve(num_samples_);
  qt_magnetizations_.reserve(num_samples_);
}

template<unsigned ISING_DIM>
template<typename T>
int Ising<ISING_DIM>::step_(int num){
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

    return 1;
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

    return 1;
  }
  //Here we only track the magnetization.
  else if constexpr(std::is_same_v<T, Wolff_type> && ISING_DIM == 3){
    int ix = disi_(rng_);
    int iy = disi_(rng_);
    int iz = disi_(rng_);
    short int target_spin = spins_[iz*N_*N_ + iy*N_ + ix]; //spin of cluster
    spins_[iz*N_*N_ + iy*N_ + ix] *= -1;
    unsigned cluster_size = 0;
    cl_wolff_queue_.push_back({ix, iy, iz});
    while(!cl_wolff_queue_.empty()){ //while we still have sites to consider
      cluster_size++;
      auto indices = cl_wolff_queue_.back();
      ix = indices[0];
      iy = indices[1];
      iz = indices[2];
      cl_wolff_queue_.pop_back();

      //look at bonds that connect to this site.
      //if bond to neighbour is occupied, flip that neighbour and add it to
      //consideration list.
      //Note: The branch ordering here has a ~10% performance impact at low temperatures.
      if(spins_[iz*N_*N_ + ((iy+1)&(N_-1))*N_ + ix] == target_spin && disd_(rng_) < cluster_growth_probability_){
        spins_[iz*N_*N_ + ((iy+1)&(N_-1))*N_ + ix] *= -1;
        cl_wolff_queue_.push_back({ix, (iy+1)&(N_-1), iz});
      }
      if(spins_[iz*N_*N_ + ((iy-1+N_)&(N_-1))*N_ + ix] == target_spin && disd_(rng_) < cluster_growth_probability_){
        spins_[iz*N_*N_ + ((iy-1+N_)&(N_-1))*N_ + ix] *= -1;
        cl_wolff_queue_.push_back({ix, (iy-1+N_)&(N_-1), iz});
      }
      if(spins_[((iz+1)&(N_-1))*N_*N_ + iy*N_ + ix] == target_spin && disd_(rng_) < cluster_growth_probability_){
        spins_[((iz+1)&(N_-1))*N_*N_ + iy*N_ + ix] *= -1;
        cl_wolff_queue_.push_back({ix, iy, (iz+1)&(N_-1)});
      }
      if(spins_[((iz-1+N_)&(N_-1))*N_*N_ + iy*N_ + ix] == target_spin && disd_(rng_) < cluster_growth_probability_){
        spins_[((iz-1+N_)&(N_-1))*N_*N_ + iy*N_ + ix] *= -1;
        cl_wolff_queue_.push_back({ix, iy, (iz-1+N_)&(N_-1)});
      }
      if(spins_[iz*N_*N_ + iy*N_ + ((ix-1+N_)&(N_-1))] == target_spin && disd_(rng_) < cluster_growth_probability_){
        spins_[iz*N_*N_ + iy*N_ + ((ix-1+N_)&(N_-1))] *= -1;
        cl_wolff_queue_.push_back({(ix-1+N_)&(N_-1), iy, iz});
      }
      if(spins_[iz*N_*N_ + iy*N_ + ((ix+1)&(N_-1))] == target_spin && disd_(rng_) < cluster_growth_probability_){
        spins_[iz*N_*N_ + iy*N_ + ((ix+1)&(N_-1))] *= -1;
        cl_wolff_queue_.push_back({(ix+1)&(N_-1), iy, iz});
      }
    }
    tracked_magnetization_ -= 2*static_cast<long int>(cluster_size)*static_cast<long int>(target_spin);
    return cluster_size;
  }
  //Here we only track the magnetization
  else if constexpr(std::is_same_v<T, SwendsenWang_type> && ISING_DIM == 3){
    //TODO: I expect this to be slow as fuck for now.
    int num_flips = 0;
    sw_reset_clusters_();
    int largest_label = 0;
    for(int iz = 0; iz < N_; ++iz){
      for(int iy = 0; iy < N_; ++iy){
        for(int ix = 0; ix < N_; ++ix){
          int index = iz*N_*N_ + iy*N_ + ix;
          int index_right = iz*N_*N_ + iy*N_ + ((ix+1)&(N_-1));
          int index_front = iz*N_*N_ + ((iy+1)&(N_-1))*N_ + ix;
          int index_above = ((iz+1)&(N_-1))*N_*N_ + iy*N_ + ix;
          //make sure we have a label
          if(sw_label_[index] == 0){
            sw_label_[index] = ++largest_label;
          }
          //check if we connect to the right
          if(spins_[index] == spins_[index_right] && disd_(rng_) < cluster_growth_probability_){
            if(sw_label_[index_right] != 0){
              sw_union_(sw_label_[index], sw_label_[index_right]);
            }
            else{
              sw_label_[index_right] = sw_label_[index];
            }
          }
          //check if we connect to the front
          if(spins_[index] == spins_[index_front] && disd_(rng_) < cluster_growth_probability_){
            if(sw_label_[index_front] != 0){
              sw_union_(sw_label_[index], sw_label_[index_front]);
            }
            else{
              sw_label_[index_front] = sw_label_[index];
            }
          }
          //check if we connect to the above
          if(spins_[index] == spins_[index_above] && disd_(rng_) < cluster_growth_probability_){
            if(sw_label_[index_above] != 0){
              sw_union_(sw_label_[index], sw_label_[index_above]);
            }
            else{
              sw_label_[index_above] = sw_label_[index];
            }
          }
        }
      }
    }
    //flip the clusters with 50% probability
    for(int i = 0; i < spins_.size(); ++i){
      //check if we have a multiplier for the corresponding cluster
      int cluster = sw_find_(sw_label_[i]);
      if(sw_flipper_[cluster] == 0){//uninitialized
        sw_flipper_[cluster] = 2 * static_cast<int>(disd_(rng_)*2.) - 1; //-1 or 1, 50% probability
      }
      if(sw_flipper_[cluster] == -1){
        ++num_flips;
        tracked_magnetization_ -= 2 * spins_[i];
      }
      spins_[i] *= sw_flipper_[cluster];
    }
    return num_flips;
  }
  else{
    throw "Ising::step_: Invalid Stepper-Dimension Combination Selected.";
  }
}

template<unsigned ISING_DIM>
template<typename T>
double Ising<ISING_DIM>::compute_energy_() const{
  if constexpr(std::is_same_v<T, MRT_type>/* && false*/){//TODO: What the hell is this?
    return tracked_energy_ / static_cast<double>(spins_.size());
  }
  else{
    double energy = 0;
    for(int iz = 0; iz < N_; ++iz){
      for(int iy = 0; iy < N_; ++iy){
        if constexpr(ISING_DIM == 2){
          energy -= J_*spins_[iz*N_+iy]*( spins_[iz*N_+((iy-1+N_)&(N_-1))] + spins_[iz*N_+((iy+1)&(N_-1))] +
                                          spins_[((iz-1+N_)&(N_-1))*N_+iy] + spins_[((iz+1)&(N_-1))*N_+iy]);
        }
        else if constexpr(ISING_DIM == 3){
          for(int ix = 0; ix < N_; ++ix){
            energy -= J_*spins_[iz*N_*N_+iy*N_+ix]*(spins_[iz*N_*N_ + iy*N_ + ((ix+1)&(N_-1))] + spins_[iz*N_*N_ + iy*N_ + ((ix-1+N_)&(N_-1))] +
                                                    spins_[iz*N_*N_ + ((iy+1)&(N_-1))*N_ + ix] + spins_[iz*N_*N_ + ((iy-1+N_)&(N_-1))*N_ + ix] +
                                                    spins_[((iz+1)&(N_-1))*N_*N_ + iy*N_ + ix] + spins_[((iz-1+N_)&(N_-1))*N_*N_ + iy*N_ + ix]);
          }
        }
      }
    }
    return energy/(2.*static_cast<double>(spins_.size()));
  }
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
  long int magnetization = 0.;
  std::for_each(spins_.begin(), spins_.end(), [&magnetization](short int s){magnetization += static_cast<long int>(s);});
  tracked_magnetization_ = magnetization;
}

template<unsigned ISING_DIM>
template<typename T>
void Ising<ISING_DIM>::fetch_sample_(){
  double energy = compute_energy_<T>();
  double magnetization = compute_magnetization_();
  energies_.push_back(energy);
  sq_energies_.push_back(energy*energy);
  magnetizations_.push_back(std::abs(magnetization));
  sq_magnetizations_.push_back(magnetization*magnetization);
  qt_magnetizations_.push_back(magnetization*magnetization*magnetization*magnetization);
}

template<unsigned ISING_DIM>
void Ising<ISING_DIM>::extract_data_(){
  double energy=0., sq_energy=0.;
  double magnetization=0., sq_magnetization=0., qt_magnetization=0.;
  //TODO: These can all use std::accumulate instead
  std::for_each(energies_.begin(), energies_.end(), [&energy](double e){energy+=e;});
  energy /= static_cast<double>(energies_.size());
  std::for_each(sq_energies_.begin(), sq_energies_.end(), [&sq_energy](double esq){sq_energy+=esq;});
  sq_energy /= static_cast<double>(sq_energies_.size());
  std::for_each(magnetizations_.begin(), magnetizations_.end(), [&magnetization](double m){magnetization+=m;});
  magnetization /= static_cast<double>(magnetizations_.size());
  std::for_each(sq_magnetizations_.begin(), sq_magnetizations_.end(), [&sq_magnetization](double msq){sq_magnetization+=msq;});
  sq_magnetization /= static_cast<double>(magnetizations_.size());
  std::for_each(qt_magnetizations_.begin(), qt_magnetizations_.end(), [&qt_magnetization](double mqt){qt_magnetization+=mqt;});
  qt_magnetization /= static_cast<double>(qt_magnetizations_.size());

  res_temperatures_.push_back(temperature_);
  res_energies_.push_back(energy);
  res_magnetizations_.push_back(magnetization);
  res_susceptibilities_.push_back((sq_magnetization-magnetization*magnetization)*static_cast<double>(spins_.size())/temperature_);
  res_capacities_.push_back((sq_energy-energy*energy)*static_cast<double>(spins_.size())/(temperature_*temperature_));
  res_binder_.push_back(1. - qt_magnetization / (3. * sq_magnetization * sq_magnetization));
}


template<unsigned ISING_DIM>
void Ising<ISING_DIM>::sw_reset_clusters_(){
  for(int i = 0; i < sw_label_.size(); ++i){
    sw_label_[i] = 0;
    sw_labels_[i] = i;
    sw_flipper_[i] = 0;
  }
}

template<unsigned ISING_DIM>
int Ising<ISING_DIM>::sw_find_(int x){
  int y = x;
  //find representative of x as y
  while(sw_labels_[y] != y){
    y = sw_labels_[y];
  }
  //chop trees
  while(sw_labels_[x] != x){
    int z = sw_labels_[x];
    sw_labels_[x] = y;
    x = z;
  }
  return y;
}

template<unsigned ISING_DIM>
void Ising<ISING_DIM>::sw_union_(const int x, const int y){
  sw_labels_[sw_find_(x)] = sw_find_(y);
}


//force compiler to produce the desired template blueprints
template void Ising<2u>::run_measurement<MRT_type>();
template void Ising<3u>::run_measurement<MRT_type>();
template void Ising<2u>::run_correlation_analysis<MRT_type>();
template void Ising<3u>::run_correlation_analysis<MRT_type>();
template void Ising<2u>::run_benchmark<MRT_type>();
template void Ising<3u>::run_benchmark<MRT_type>();
template void Ising<3u>::run_measurement<Wolff_type>();
template void Ising<3u>::run_correlation_analysis<Wolff_type>();
template void Ising<3u>::run_benchmark<Wolff_type>();
template void Ising<3u>::run_measurement<SwendsenWang_type>();
template void Ising<3u>::run_correlation_analysis<SwendsenWang_type>();
template void Ising<3u>::run_benchmark<SwendsenWang_type>();
template class Ising<2u>;
template class Ising<3u>;
