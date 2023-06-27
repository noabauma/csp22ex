/* sanity_check.cpp 2021
 * Author: Pascal Engeler
 * engelerp@phys.ethz.ch
 */

#include "Ising.hpp"
#include <iostream>
#include <string>
#include <array>
#include <cmath>

int main(){
  std::cout << "Running sanity check" << std::endl;
  std::vector<double> temperatures;
  //we cluster the target temperatures around the critical temperature
  for(double t = -1.5; t < 1.5; t += 0.1){
    temperatures.push_back(4.51 + 3.*(t)*(t)*(t)/5.3333);
  }
  const unsigned SIZE = 8;
  const unsigned numspins = SIZE*SIZE*SIZE;

  //number of Wolff steps needed to get about one sweep
  //ordered like the temperatures
  std::vector<unsigned> updates_per_sweep =
  { 2,2,2,2,2,2,2,3,3,4,5,6,6,6,6,6,6,6,7,7,8,10,13,17,24,32,43,57,72,89 };
  std::vector<unsigned> updates_per_sweep_sw =
  { 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2 };


  Ising<3> sys3D_Wolff (SIZE, 1.);
  Ising<3> sys3D_SW (SIZE, 1.);
  Ising<3> sys3D_MRT (SIZE, 1.);

  unsigned thermalize_sweeps = 60;
  unsigned num_samples = 20000;
  unsigned sample_stride = 6;

  std::cout << "Generating Wolff, Swendsen-Wang and MRT sanity check data at 8x8x8..." << std::endl;
  for(int i = 0; i < temperatures.size(); ++i){
    std::cout << temperatures[i] << std::flush << " ";
    sys3D_Wolff.set_measurement_parameters(temperatures[i], thermalize_sweeps*updates_per_sweep[i], num_samples, sample_stride*updates_per_sweep[i]);
    sys3D_Wolff.run_measurement<Wolff_type>();
    sys3D_SW.set_measurement_parameters(temperatures[i], thermalize_sweeps*updates_per_sweep_sw[i], num_samples, sample_stride*updates_per_sweep_sw[i]);
    sys3D_SW.run_measurement<SwendsenWang_type>();
    sys3D_MRT.set_measurement_parameters(temperatures[i], thermalize_sweeps*numspins, num_samples, sample_stride*numspins);
    sys3D_MRT.run_measurement<MRT_type>();
  }
  std::cout << std::endl;
  sys3D_Wolff.save_data("sanity_check_wolff.csv");
  sys3D_SW.save_data("sanity_check_sw.csv");
  sys3D_MRT.save_data("sanity_check_mrt.csv");


  //More Binder cumulant data for Wolff
  std::cout << "Generating Wolff Binder Cumulant data..." << std::endl;
  std::array<unsigned, 2> sizes {4,16};
  for(auto size: sizes){
    std::cout << "L = " << size << std::endl;
    Ising<3> sys3D (size, 1.);
    for(int i = 0; i < temperatures.size(); ++i){
      std::cout << temperatures[i] << std::flush << " ";
      sys3D.set_measurement_parameters(temperatures[i], thermalize_sweeps*updates_per_sweep[i], num_samples, sample_stride*updates_per_sweep[i]);
      sys3D.run_measurement<Wolff_type>();
    }
    std::cout << std::endl;
    sys3D.save_data("sanity_check_wolff_bc_"+std::to_string(size)+".csv");
  }

  std::cout << std::endl << "Sanity check complete." << std::endl << std::endl;

  return 0;
}
