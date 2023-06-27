/* correlation_analysis.cpp 2021
 * Author: Pascal Engeler
 * engelerp@phys.ethz.ch
 */

#include "Ising.hpp"
#include <iostream>
#include <string>

int main(){
  const std::vector<unsigned> sizes = {8, 16, 32}; //system sizes we want to consider
  std::vector<double> temperatures = {4.,5.,6}; //temperatures we want to consider
  for(auto SZE: sizes){
    std::cout << "Size: " << SZE << std::endl;
    Ising<3> sys3D (SZE, 1.); //construct ising
    sys3D.seed(42); //seed the rng
    std::cout << "T = ";
    for(auto t: temperatures){
      std::cout << t << "\t" << std::flush;
      sys3D.set_measurement_parameters(t, 20*SZE*SZE*SZE, 20*SZE*SZE*SZE, 1); //setup measurement details
      sys3D.run_correlation_analysis<MRT_type>(); //run correlation analysis
    }
    std::cout << std::endl << "Saving..." << std::flush;
    std::string sstr = std::to_string(SZE);
    sys3D.save_correlation_data("corrdata_" + sstr + "x" + sstr + "x" + sstr +  ".csv"); //save the data
    std::cout << std::endl;
  }

  return 0;
}
