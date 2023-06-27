/* speed.cpp 2021
 * Author: Pascal Engeler
 * engelerp@phys.ethz.ch
 */

#include "Ising.hpp"
#include <iostream>
#include <string>
#include <array>
#include <cmath>

int main(){
  std::cout << "Running monte carlo speed calculation" << std::endl;
  std::vector<double> temperatures = {3.51, 4.51, 5.51};
  const unsigned SIZE = 16;
  unsigned numspins = SIZE*SIZE*SIZE;

  //number of Wolff steps needed to get about one sweep
  //ordered like the temperatures
  std::vector<unsigned> updates_per_sweep = {2, 10, 500};
  std::vector<unsigned> updates_per_sweep_sw = {2, 2, 2};


  //Wolff for 3 different temperatures
  std::cout << "Wolff for 3 different temperatures..." << std::endl;
  Ising<3> sys3D_Wolff (SIZE, 1.);
  for(int i = 0; i < temperatures.size(); ++i){
    std::cout << temperatures[i] << std::flush << " ";
    sys3D_Wolff.set_measurement_parameters(temperatures[i], 100*updates_per_sweep[i], 20000*updates_per_sweep[i], 0);
    sys3D_Wolff.run_benchmark<Wolff_type>();
    sys3D_Wolff.set_measurement_parameters(temperatures[i], 100*updates_per_sweep[i], 20000, 0);
    sys3D_Wolff.run_correlation_analysis<Wolff_type>();
  }
  sys3D_Wolff.save_correlation_data("correlations_wolff_16.csv");
  sys3D_Wolff.save_benchmark_data("benchmark_wolff_16.csv");

  //SW for 3 different temperatures
  std::cout << "SW for 3 different temperatures..." << std::endl;
  Ising<3> sys3D_SW (SIZE, 1.);
  for(int i = 0; i < temperatures.size(); ++i){
    std::cout << temperatures[i] << std::flush << " ";
    sys3D_SW.set_measurement_parameters(temperatures[i], 100*updates_per_sweep_sw[i], 20000*updates_per_sweep_sw[i], 0);
    sys3D_SW.run_benchmark<SwendsenWang_type>();
    sys3D_SW.set_measurement_parameters(temperatures[i], 100*updates_per_sweep_sw[i], 20000, 0);
    sys3D_SW.run_correlation_analysis<SwendsenWang_type>();
  }
  sys3D_SW.save_correlation_data("correlations_sw_16.csv");
  sys3D_SW.save_benchmark_data("benchmark_sw_16.csv");

  //MRT for 3 different temperatures
  std::cout << "MRT for 3 different temperatures..." << std::endl;
  Ising<3> sys3D_MRT (SIZE, 1.);
  for(int i = 0; i < temperatures.size(); ++i){
    std::cout << temperatures[i] << std::flush << " ";
    sys3D_MRT.set_measurement_parameters(temperatures[i], 100*numspins, 20000*numspins, 0);
    sys3D_MRT.run_benchmark<MRT_type>();
    sys3D_MRT.set_measurement_parameters(temperatures[i], 100*numspins, 20000, 0);
    sys3D_MRT.run_correlation_analysis<MRT_type>();
  }
  sys3D_MRT.save_correlation_data("correlations_mrt_16.csv");
  sys3D_MRT.save_benchmark_data("benchmark_mrt_16.csv");

  //Single temperature, several sizes
  std::cout << "Several system sizes at Tc..." << std::endl;
  std::vector<unsigned> sizes = {2, 4, 8, 16, 32};
  for(auto size: sizes){
    std::cout << "L = " << size << std::endl;
    std::cout << " Wolff "<< std::flush;
    numspins = size*size*size;
    Ising<3> sys3D_Wolff (size, 1.);
    sys3D_Wolff.set_measurement_parameters(temperatures[1], 100*updates_per_sweep[1], 20000*updates_per_sweep[1], 0);
    sys3D_Wolff.run_benchmark<Wolff_type>();
    sys3D_Wolff.set_measurement_parameters(temperatures[1], 100*updates_per_sweep[1], 20000, 0);
    sys3D_Wolff.run_correlation_analysis<Wolff_type>();
    std::string sstr = std::to_string(size);
    sys3D_Wolff.save_correlation_data("tc_correlations_wolff_"+sstr+".csv");
    sys3D_Wolff.save_benchmark_data("tc_benchmark_wolff_"+sstr+".csv");
    std::cout << "\n SW " << std::flush;
    Ising<3> sys3D_SW (size, 1.);
    sys3D_SW.set_measurement_parameters(temperatures[1], 100*updates_per_sweep_sw[1], 20000*updates_per_sweep_sw[1], 0);
    sys3D_SW.run_benchmark<SwendsenWang_type>();
    sys3D_SW.set_measurement_parameters(temperatures[1], 100*updates_per_sweep_sw[1], 20000, 0);
    sys3D_SW.run_correlation_analysis<SwendsenWang_type>();
    sys3D_SW.save_correlation_data("tc_correlations_sw_"+sstr+".csv");
    sys3D_SW.save_benchmark_data("tc_benchmark_sw_"+sstr+".csv");
    std::cout << "\n MRT " << std::flush;
    Ising<3> sys3D_MRT (size, 1.);
    sys3D_MRT.set_measurement_parameters(temperatures[1], 100*numspins, 20000*numspins, 0);
    sys3D_MRT.run_benchmark<MRT_type>();
    sys3D_MRT.set_measurement_parameters(temperatures[1], 100*numspins, 20000, 0);
    sys3D_MRT.run_correlation_analysis<MRT_type>();
    sys3D_MRT.save_correlation_data("tc_correlations_mrt_"+sstr+".csv");
    sys3D_MRT.save_benchmark_data("tc_benchmark_mrt_"+sstr+".csv");
  }


  std::cout << std::endl << "Monte carlo speed calculation finished." << std::endl;
  return 0;
}
