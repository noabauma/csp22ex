/* main.cpp 2021
 * Author: Pascal Engeler
 * engelerp@phys.ethz.ch
 */

#include "Ising.hpp"
#include <iostream>
#include <string>
#include <array>
#include <cmath>

int main(){
  const std::vector<unsigned> sizes = {2, 4, 8, 16, 32};
  std::vector<double> temperatures;
  //we cluster the target temperatures around the critical temperature
  for(double t = -1.5; t < 1.5; t += 0.1){
    temperatures.push_back(4.51 + 3.*(t)*(t)*(t)/5.3333);
  }
  unsigned num_samples = 4000; //we take 4000 samples per temperature
  unsigned thermalize_sweeps = 30; //number of sweeps we use to thermalize
  std::vector<unsigned> sample_strides; //number of sweeps we drop between successive samples
  for(auto t: temperatures){
    if(std::abs(t - 4.51) > 0.3){ //outside of the critical region
      sample_strides.push_back(5);
    }
    else{ //within the critical region
      sample_strides.push_back(12); //we just use the largest value for all of them, so we don't need to setup individually
    }
  }

  for(auto SIZE: sizes){
    std::cout << "Linear size: " << SIZE << std::endl;
    Ising<3> sys3D (SIZE, 1.);
    sys3D.seed(42);
    for(int i = 0; i < temperatures.size(); ++i){
      std::cout << "\t" << temperatures[i] << std::flush;
      sys3D.set_measurement_parameters(temperatures[i], thermalize_sweeps*SIZE*SIZE*SIZE, num_samples, sample_strides[i]*SIZE*SIZE*SIZE);
      sys3D.run_measurement<MRT_type>();
    }
    std::cout << std::endl;
    sys3D.save_data("data_" + std::to_string(3)+"D_"+std::to_string(SIZE)+".csv");
  }

  return 0;
}
