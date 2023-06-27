/* main.cpp 2021
 * Author: Pascal Engeler
 * engelerp@phys.ethz.ch
 */

#include "Ising.hpp"
#include <vector>
#include <fstream>
#include <string>
#include <iostream>

int main(){
  const std::vector<double> train_temperatures = {
              0.        , 0.22222222, 0.44444444, 0.66666667, 0.88888889,
              1.11111111, 1.33333333, 1.55555556, 1.77777778, 2.,
              2.6       , 2.81111111, 3.02222222, 3.23333333, 3.44444444,
              3.65555556, 3.86666667, 4.07777778, 4.28888889, 4.5       };

  const std::vector<double> verify_temperatures = {
              0.        , 0.04545455, 0.09090909, 0.13636364, 0.18181818,
              0.22727273, 0.27272727, 0.31818182, 0.36363636, 0.40909091,
              0.45454545, 0.5       , 0.54545455, 0.59090909, 0.63636364,
              0.68181818, 0.72727273, 0.77272727, 0.81818182, 0.86363636,
              0.90909091, 0.95454545, 1.        , 1.04545455, 1.09090909,
              1.13636364, 1.18181818, 1.22727273, 1.27272727, 1.31818182,
              1.36363636, 1.40909091, 1.45454545, 1.5       , 1.54545455,
              1.59090909, 1.63636364, 1.68181818, 1.72727273, 1.77272727,
              1.81818182, 1.86363636, 1.90909091, 1.95454545, 2.        ,
              2.04545455, 2.09090909, 2.13636364, 2.18181818, 2.22727273,
              2.27272727, 2.31818182, 2.36363636, 2.40909091, 2.45454545,
              2.5       , 2.54545455, 2.59090909, 2.63636364, 2.68181818,
              2.72727273, 2.77272727, 2.81818182, 2.86363636, 2.90909091,
              2.95454545, 3.        , 3.04545455, 3.09090909, 3.13636364,
              3.18181818, 3.22727273, 3.27272727, 3.31818182, 3.36363636,
              3.40909091, 3.45454545, 3.5       , 3.54545455, 3.59090909,
              3.63636364, 3.68181818, 3.72727273, 3.77272727, 3.81818182,
              3.86363636, 3.90909091, 3.95454545, 4.        , 4.04545455,
              4.09090909, 4.13636364, 4.18181818, 4.22727273, 4.27272727,
              4.31818182, 4.36363636, 4.40909091, 4.45454545, 4.5       };

  using stepper = MRT_type;
  const size_t N = 32;
  const size_t samp_per_temp = 1000;
  const size_t thermalize_steps = N*N*50;
  const size_t decorrelat_steps = N*N*10;
  const std::string train_temps_file = "temps_train.txt";
  const std::string train_spins_file = "spins_train.txt";
  const std::string verify_temps_file = "temps_verify.txt";
  const std::string verify_spins_file = "spins_verify.txt";
  std::ofstream file;

  Ising<2> ising_model_train (N,1.);
  std::vector<short int> spinvec (N * N * samp_per_temp, 0);

  //Generate training data
  std::cout << "Generating Training Data\nT: ";
  for(auto T: train_temperatures){
    std::cout << T << " " << std::flush;
    ising_model_train.set_measurement_parameters(T, 0, 0, 0);
    //thermalize
    ising_model_train.step_N<stepper>(thermalize_steps);
    //sample
    for(size_t i = 0; i < samp_per_temp; ++i){
      //decorrelate
      ising_model_train.step_N<stepper>(decorrelat_steps);
      //fetch sample and copy it over
      const std::vector<short int>& spins = ising_model_train.get_configuration();
      std::copy(spins.begin(), spins.end(), spinvec.begin()+N*N*i);
    }
    //write samples for that temperature to file
    file.open(train_spins_file, std::ios::app);
    if(!file){
      std::cout << "Can't open file "+train_spins_file;
      return -1;
    }
    for(size_t i = 0; i < spinvec.size(); ++i){
      if(((i % (N*N)) == 0) && i != 0)
        file << "\n";
      file << spinvec[i] << " ";
    }
    file << "\n";
    file.close();
    //write temperature to file
    file.open(train_temps_file, std::ios::app);
    if(!file){
      std::cout << "Can't open file "+train_temps_file;
      return -1;
    }
    for(size_t i = 0; i < samp_per_temp; ++i){
      file << T << "\n";
    }
    file.close();
  }

  Ising<2> ising_model_predict (N,1.);
  //Generate verification data
  std::cout << "\nGenerating Verification Data\nT: ";
  for(auto T: verify_temperatures){
    std::cout << T << " " << std::flush;
    ising_model_predict.set_measurement_parameters(T, 0, 0, 0);
    //thermalize
    ising_model_predict.step_N<stepper>(thermalize_steps);
    //sample
    for(size_t i = 0; i < samp_per_temp; ++i){
      //decorrelate
      ising_model_predict.step_N<stepper>(decorrelat_steps);
      //fetch sample and copy it over
      const std::vector<short int>& spins = ising_model_predict.get_configuration();
      std::copy(spins.begin(), spins.end(), spinvec.begin()+N*N*i);
    }
    //write samples for that temperature to file
    file.open(verify_spins_file, std::ios::app);
    if(!file){
      std::cout << "Can't open file "+verify_spins_file;
      return -1;
    }
    for(size_t i = 0; i < spinvec.size(); ++i){
      if(((i % (N*N)) == 0) && i != 0)
        file << "\n";
      file << spinvec[i] << " ";
    }
    file << "\n";
    file.close();
    //write temperature to file
    file.open(verify_temps_file, std::ios::app);
    if(!file){
      std::cout << "Can't open file "+verify_temps_file;
      return -1;
    }
    for(size_t i = 0; i < samp_per_temp; ++i){
      file << T << "\n";
    }
    file.close();
  }

  return 0;
}
