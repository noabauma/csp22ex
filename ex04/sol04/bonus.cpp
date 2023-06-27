#include "Ising.hpp"
#include <iostream>
#include <string>
#include <array>
#include <cmath>

int main(){
  std::vector<double> temperatures;
  //we cluster the target temperatures around the critical temperature
  for(double t = -1.5; t < 1.5; t += 0.1){
    temperatures.push_back(4.51 + 3.*(t)*(t)*(t)/5.3333);
  }
  const unsigned SIZE = 16;
  unsigned numspins = SIZE*SIZE*SIZE;

  //number of Wolff steps needed to get about one sweep
  //ordered like the temperatures
  std::vector<unsigned> updates_per_sweep =
    {2,2,2,2,2,2,2,3,3,4,6,7,10,11,11,12,12,13,16,16,26,37,83,166,193,249,400,468,586,721};
  std::vector<unsigned> updates_per_sweep_sw =
    {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};

  //Wolff for different temperatures
  Ising<3> sys3D_Wolff (SIZE, 1.);
  for(int i = 0; i < temperatures.size(); ++i){
    std::cout << temperatures[i] << std::flush << " ";
    sys3D_Wolff.set_measurement_parameters(temperatures[i], 300*updates_per_sweep[i], 50000*updates_per_sweep[i], 0);
    sys3D_Wolff.run_benchmark<Wolff_type>();
    sys3D_Wolff.set_measurement_parameters(temperatures[i], 300*updates_per_sweep[i], 50000, 0);
    sys3D_Wolff.run_correlation_analysis<Wolff_type>();
  }

  sys3D_Wolff.save_correlation_data("bonus_correlations_wolff_16.csv");
  sys3D_Wolff.save_benchmark_data("bonus_benchmark_wolff_16.csv");

  //SW for different temperatures
  Ising<3> sys3D_SW (SIZE, 1.);
  for(int i = 0; i < temperatures.size(); ++i){
    std::cout << temperatures[i] << std::flush << " ";
    sys3D_SW.set_measurement_parameters(temperatures[i], 300*updates_per_sweep_sw[i], 50000*updates_per_sweep_sw[i], 0);
    sys3D_SW.run_benchmark<SwendsenWang_type>();
    sys3D_SW.set_measurement_parameters(temperatures[i], 300*updates_per_sweep_sw[i], 50000, 0);
    sys3D_SW.run_correlation_analysis<SwendsenWang_type>();
  }

  sys3D_SW.save_correlation_data("bonus_correlations_sw_16.csv");
  sys3D_SW.save_benchmark_data("bonus_benchmark_sw_16.csv");

  //MRT for different temperatures
  Ising<3> sys3D_MRT (SIZE, 1.);
  for(int i = 0; i < temperatures.size(); ++i){
    std::cout << temperatures[i] << std::flush << " ";
    sys3D_MRT.set_measurement_parameters(temperatures[i], 300*numspins, 50000*numspins, 0);
    sys3D_MRT.run_benchmark<MRT_type>();
    sys3D_MRT.set_measurement_parameters(temperatures[i], 300*numspins, 50000, 0);
    sys3D_MRT.run_correlation_analysis<MRT_type>();
  }
  sys3D_MRT.save_correlation_data("bonus_correlations_mrt_16.csv");
  sys3D_MRT.save_benchmark_data("bonus_benchmark_mrt_16.csv");


  return 0;
}
