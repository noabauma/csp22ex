//#include <md_box_nh_naive.hpp>
#include <md_box_nh.hpp>
#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>

//Struct to hold system data
struct Results {
  Results(double energy, double temperature, double xi, double ekin, double maxvel): energy(energy), temperature(temperature), xi(xi), ekin(ekin), maxvel(maxvel) {}

  double energy;
  double temperature;
  double xi;
  double ekin;
  double maxvel;
};

//Struct to hold speed data
template <size_t N>
struct SpeedData {
  SpeedData(size_t timestep, std::array<double, N> speeds): timestep(timestep), speeds(speeds) {}

  size_t timestep;
  std::array<double, N> speeds;
};

std::ostream& operator<<(std::ostream& os, const Results& res) {
  os << res.energy << " " << res.temperature << " " << res.xi << " " << res.ekin << " " << res.maxvel;
  return os;
}

template <size_t N>
std::ostream& operator<<(std::ostream& os, const SpeedData<N>& sd) {
  os << sd.timestep;
  for (auto speed : sd.speeds) {
    os << " " << speed;
  }
  return os;
}

//Helper to generate filenames
std::string double_to_string(double x) {
  std::string inter = std::to_string(x);
  std::string before = "";
  std::string after = "";
  bool comma_seen = false;
  //split before and after decimal point
  for (size_t i = 0; i < inter.size(); ++i) {
    if (inter[i] == '.') {
      comma_seen = true;
      continue;
    }
    if (!comma_seen) {
      before += inter[i];
    }
    else {
      after += inter[i];
    }
  }
  if (after.size() > 0) {
    //erase trailing zeros
    int current_i = after.size() - 1;
    while (current_i >= 0 && after[current_i] == '0') {
      --current_i;
    }
    after.resize(current_i + 1);
  }
  
  if (after.size() > 0) {
    return before + "d" + after;
  }
  else {
    return before;
  }
}

//You may want to specify an abnormally large stack size for this one
int main(int argc, char* argv[]){
  if(argc != 3){
    std::cout << "Usage: ./" << argv[0] << " Q num_steps" << std::endl;
    return 1;
  }
  double Q = std::atof(argv[1]);
  int num_steps = std::atoi(argv[2]);
  if(num_steps < 0){
    std::cout << "Invalid number of steps: " << num_steps << std::endl;
    return 1;
  }

  const size_t N = 125;
  MD_Box<double, N, 3> box (0.001, 41, Q, 2.);

  std::vector<Results> resvec;
  resvec.reserve(num_steps);
  std::vector<SpeedData<N>> sdvec;
  sdvec.reserve(10);
  unsigned sbin_stride = num_steps / 10;

  for (size_t i = 0; i < num_steps; ++i) {
    double energy = box.step_all();
    resvec.push_back({ energy, box.temperature(), box.xi(), box.Ekin(), box.max_velocity() });
    if (i % sbin_stride == 0) {
      sdvec.push_back(SpeedData<N>(i, box.speeds()));
    }
  }

  //write system data
  std::string filename = "data_Q" + double_to_string(Q) + ".txt";
  std::ofstream file(filename);

  try {
    for (auto r : resvec) {
      file << r << "\n";
    }
    file.flush();
    file.close();
  }
  catch (...) {
    std::cout << "Error writing to file " << filename << std::endl;
    return 1;
  }

  //write speed data
  filename = "data_Q" + double_to_string(Q) + "_speeds.txt";
  std::ofstream file2(filename);

  try {
    for (auto r : sdvec) {
      file2 << r << "\n";
    }
    file2.flush();
    file2.close();
  }
  catch (...) {
    std::cout << "Error writing to file " << filename << std::endl;
    return 1;
  }

  return 0;
}
