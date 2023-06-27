/* Computational Statistical Physics 2021
 * Exercise 8
 * Author: Pascal Engeler, engelerp@phys.ethz.ch
 */
#ifndef MD_BOX_HPP_INCLUDED
#define MD_BOX_HPP_INCLUDED
#include <vec3.hpp>
#include <array>
#include <lennard_jones.hpp>
#include <algorithm>
#include <random>
#include <iostream>

//Note that here, L is in units of 2.5 * sigma, i.e. it's the number of cells per dimension
//Cells have linear size 2.5, equal to the cutoff, and sigma=1.
template <typename real_t, size_t N, size_t L>
class MD_Box{
public:

  MD_Box(real_t dt = 0.01, unsigned seed = 42, real_t Q = 100, real_t T = 100) : potential_(1,1), force_(1,1), L_(2.5 * L), NM_(L), dt_(dt), xi_0_(0), Q_(Q), temperature_(T)
  {
    std::mt19937_64 rng (seed);
    std::normal_distribution<real_t> dis_normal (real_t(0), real_t(1));
    //fill in previous positions by choosing randomly from well spaced positions
    std::vector<Vec3<real_t> > possible_positions;
    for (real_t x = 0; x < L_; x += 1.5) {
      for (real_t y = 0; y < L_; y += 1.5) {
        for (real_t z = 0; z < L_; z += 1.5) {
          possible_positions.push_back({ x,y,z });
        }
      }
    }
    if (N > possible_positions.size()) {
      std::cout << "Can't place that many particles in this system. Maximum: " << possible_positions.size() << std::endl;
      exit(1);
    }
    std::shuffle(possible_positions.begin(), possible_positions.end(), rng);
    for (size_t i = 0; i < particle_positions_0_.size(); ++i) {
      particle_positions_0_[i] = possible_positions[i];
      periodize_position_(particle_positions_0_[i]); //this should be 0x90
    }

    //generate initial velocities
    for(size_t i = 0; i < particle_velocities_0_.size(); ++i){
      Vec3<real_t> velocity(dis_normal(rng), dis_normal(rng), dis_normal(rng));
      particle_velocities_0_[i] = velocity;
      //particle_velocities_0_[i] = {0, 0, 0};
    }

    //reset FIRST_
    for (size_t i = 0; i < FIRST_.size(); ++i) {
      FIRST_[i] = -1;
    }
    //assign the particles into cells. Note that cell (i + NM_ * j + NM_ * NM_ * k) has its lower left corner at (2.5 * i, 2.5 * j, 2.5 * k).
    for (size_t p_i = 0; p_i < N; ++p_i) {
      add_particle_to_cell_(p_i, ijk_to_cell_index_(pos_to_ijk_(particle_positions_0_[p_i])));
    }
  }



  //We iterate over all cells twice:
  //0 - update kinetic energy 0
  //1 - substep 1 for all cells
  //update kinetic energy
  real_t step_all() noexcept {
    std::vector<size_t> storage;
    storage.reserve(800);
    for (size_t k = 0; k < NM_; ++k) {
      for (size_t j = 0; j < NM_; ++j) {
        for (size_t i = 0; i < NM_; ++i) {
          if (FIRST_[ijk_to_cell_index_({ i,j,k })] == -1) {
            continue;
          }
          storage.clear();
          substep_cell_1_({ i,j,k }, storage);
        }
      }
    }
    compute_Ekin_0_();
    compute_xi_1_();
    compute_Ekin_1_();
    compute_xi_2_();
    //swap new particle positions into place
    std::swap(particle_positions_0_, particle_positions_1_);
    //and re-periodize them
    for (size_t i = 0; i < N; ++i) {
      periodize_position_(particle_positions_0_[i]);
    }
    //fix particle-to-cell assignment
    reassign_particles_();
    //now do substep 2
    for (size_t k = 0; k < NM_; ++k) {
      for (size_t j = 0; j < NM_; ++j) {
        for (size_t i = 0; i < NM_; ++i) {
          if (FIRST_[ijk_to_cell_index_({ i,j,k })] == -1) {
            continue;
          }
          storage.clear();
          substep_cell_2_({ i,j,k }, storage);
        }
      }
    }
    //swap rest of new stuff into place
    std::swap(particle_velocities_0_, particle_velocities_2_);
    xi_0_ = xi_2_;
    Ekin_0_ = Ekin_1_;

    return energy();
  }

  const std::array<Vec3<real_t>, N>& positions() {
    return particle_positions_0_;
  }

  real_t max_velocity() {
    real_t maxvsq = 0;
    for(auto v: particle_velocities_0_){
      if(v.rsq() > maxvsq){
        maxvsq = v.rsq();
      }
    }
    return std::sqrt(maxvsq);
  }

  std::array<double, N> speeds() {
    std::array<double, N> spds;
    std::transform(particle_velocities_0_.begin(), particle_velocities_0_.end(), spds.begin(), [](auto v) { return v.r(); });
    return spds;
  }

  //returns system energy at positions 0 and velocities 0
  real_t energy(){
    real_t energy_pot(0), energy_kin(0);
    std::vector<size_t> storage;
    storage.reserve(800);
    for (size_t k = 0; k < NM_; ++k) {
      for (size_t j = 0; j < NM_; ++j) {
        for (size_t i = 0; i < NM_; ++i) {
          if (FIRST_[ijk_to_cell_index_({ i,j,k })] == -1) {
            continue;
          }
          storage.clear();
          const Vec3<size_t> ijk = { i, j, k };
          int particle = FIRST_[ijk_to_cell_index_(ijk)];
          get_all_in_cells_(neighbour_cell_indices_(ijk), storage);
          while(particle != -1){
            Vec3<real_t> current_particle = particle_positions_0_[particle];
            for(auto p: storage){//calculate force on particle
              if(p == particle){
                continue;
              }
              else{
                Vec3<real_t> distance = particle_positions_0_[p] - current_particle;
                periodize_distance_(distance);
                energy_pot += potential_(distance.r());
              }
            }
            particle = LIST_[particle];
          }
        }
      }
    }

    //do kinetic energy here
    for(auto v: particle_velocities_0_){
      energy_kin += v.rsq();
    }
    Ekin_0_ = 0.5 * energy_kin;

    return 0.5 * (energy_kin + energy_pot);
  }

  //This assumes the correct kinetic energy is in Ekin_0_
  real_t temperature(){
    return 2. / (3. * (N - 1.)) * Ekin_0_;
  }

  real_t dt() { return dt_; }

  real_t xi() { return xi_0_; }

  real_t Ekin() { return Ekin_0_; }

private:
  inline void periodize_distance_(Vec3<real_t>& distance){ //This can potentially be more efficient without branching
    distance.x() = (std::abs(distance.x()) < (0.5 * L_)) ? distance.x() : std::copysign(L_ - std::abs(distance.x()), -distance.x());
    distance.y() = (std::abs(distance.y()) < (0.5 * L_)) ? distance.y() : std::copysign(L_ - std::abs(distance.y()), -distance.y());
    distance.z() = (std::abs(distance.z()) < (0.5 * L_)) ? distance.z() : std::copysign(L_ - std::abs(distance.z()), -distance.z());
  }
  inline void periodize_position_(Vec3<real_t>& pos){
    if (pos.x() >= L_) {
      pos.x() -= L_;
    }
    else if (pos.x() < 0) {
      pos.x() += L_;
    }
    if (pos.y() >= L_) {
      pos.y() -= L_;
    }
    else if (pos.y() < 0) {
      pos.y() += L_;
    }
    if (pos.z() >= L_) {
      pos.z() -= L_;
    }
    else if (pos.z() < 0) {
      pos.z() += L_;
    }
  }
  inline Vec3<size_t> pos_to_ijk_(const Vec3<real_t>& p) {
    Vec3<size_t> ret_vec = { static_cast<size_t>(p.x() / real_t(2.5)), static_cast<size_t>(p.y() / real_t(2.5)), static_cast<size_t>(p.z() / real_t(2.5)) };
    return ret_vec;
  }
  inline size_t ijk_to_cell_index_(const Vec3<size_t>& ijk) {
    return ijk.x() + ijk.y() * NM_ + ijk.z() * NM_ * NM_;
  }
  //Modulo operation can potentially be made more efficient. alternatives: branching, bitwise for powers of two.
  //profile these, see what's what.
  inline std::array<Vec3<size_t>, 27> neighbour_cell_indices_(const Vec3<size_t>& ijk) {
    int i = ijk.x();
    int j = ijk.y();
    int k = ijk.z();
    std::array<Vec3<size_t>, 27> ret_array = { { //std::array requires aggregate initialization of array inside (guh)
      {static_cast<size_t>((i - 1 + L) % L), static_cast<size_t>((j - 1 + L) % L),  static_cast<size_t>((k - 1 + L) % L)},
      {static_cast<size_t>(i),           static_cast<size_t>((j - 1 + L) % L),  static_cast<size_t>((k - 1 + L) % L)},
      {static_cast<size_t>((i + 1) % L), static_cast<size_t>((j - 1 + L) % L),  static_cast<size_t>((k - 1 + L) % L)},
      {static_cast<size_t>((i - 1 + L) % L), static_cast<size_t>(j),            static_cast<size_t>((k - 1 + L) % L)},
      {static_cast<size_t>(i),           static_cast<size_t>(j),            static_cast<size_t>((k - 1 + L) % L)},
      {static_cast<size_t>((i + 1) % L), static_cast<size_t>(j),            static_cast<size_t>((k - 1 + L) % L)},
      {static_cast<size_t>((i - 1 + L) % L), static_cast<size_t>((j + 1) % L),  static_cast<size_t>((k - 1 + L) % L)},
      {static_cast<size_t>(i),           static_cast<size_t>((j + 1) % L),  static_cast<size_t>((k - 1 + L) % L)},
      {static_cast<size_t>((i + 1) % L), static_cast<size_t>((j + 1) % L),  static_cast<size_t>((k - 1 + L) % L)},
      {static_cast<size_t>((i - 1 + L) % L), static_cast<size_t>((j - 1 + L) % L),  static_cast<size_t>(k)},
      {static_cast<size_t>(i),           static_cast<size_t>((j - 1 + L) % L),  static_cast<size_t>(k)},
      {static_cast<size_t>((i + 1) % L), static_cast<size_t>((j - 1 + L) % L),  static_cast<size_t>(k)},
      {static_cast<size_t>((i - 1 + L) % L), static_cast<size_t>(j),            static_cast<size_t>(k)},
      {static_cast<size_t>(i),           static_cast<size_t>(j),            static_cast<size_t>(k)},
      {static_cast<size_t>((i + 1) % L), static_cast<size_t>(j),            static_cast<size_t>(k)},
      {static_cast<size_t>((i - 1 + L) % L), static_cast<size_t>((j + 1) % L),  static_cast<size_t>(k)},
      {static_cast<size_t>(i),           static_cast<size_t>((j + 1) % L),  static_cast<size_t>(k)},
      {static_cast<size_t>((i + 1) % L), static_cast<size_t>((j + 1) % L),  static_cast<size_t>(k)},
      {static_cast<size_t>((i - 1 + L) % L), static_cast<size_t>((j - 1 + L) % L),  static_cast<size_t>((k + 1) % L)},
      {static_cast<size_t>(i),           static_cast<size_t>((j - 1 + L) % L),  static_cast<size_t>((k + 1) % L)},
      {static_cast<size_t>((i + 1) % L), static_cast<size_t>((j - 1 + L) % L),  static_cast<size_t>((k + 1) % L)},
      {static_cast<size_t>((i - 1 + L) % L), static_cast<size_t>(j),            static_cast<size_t>((k + 1) % L)},
      {static_cast<size_t>(i),           static_cast<size_t>(j),            static_cast<size_t>((k + 1) % L)},
      {static_cast<size_t>((i + 1) % L), static_cast<size_t>(j),            static_cast<size_t>((k + 1) % L)},
      {static_cast<size_t>((i - 1 + L) % L), static_cast<size_t>((j + 1) % L),  static_cast<size_t>((k + 1) % L)},
      {static_cast<size_t>(i),           static_cast<size_t>((j + 1) % L),  static_cast<size_t>((k + 1) % L)},
      {static_cast<size_t>((i + 1) % L), static_cast<size_t>((j + 1) % L),  static_cast<size_t>((k + 1) % L)}} };
    return ret_array;
  }

  inline bool remove_particle_from_cell_(size_t p_i, size_t c_i) {
    if (FIRST_[c_i] == p_i) {
      FIRST_[c_i] = LIST_[p_i];
      cell_of_particle_[p_i] = -1;
      return true;
    }
    int val = FIRST_[c_i];
    int pre = -1;
    while (val != p_i) {
      if (val == -1) [[unlikely]] {
        return false;
      }
      pre = val;
      val = LIST_[val];
    }
    LIST_[pre] = LIST_[p_i];
    cell_of_particle_[p_i] = -1;
    return true;
  }
  inline void add_particle_to_cell_(size_t p_i, size_t c_i) {
    if (FIRST_[c_i] == -1) {
      FIRST_[c_i] = p_i;
      LIST_[p_i] = -1;
      cell_of_particle_[p_i] = c_i;
      return;
    }
    else {
      LIST_[p_i] = FIRST_[c_i];
      FIRST_[c_i] = p_i;
      cell_of_particle_[p_i] = c_i;
      return;
    }
  }
  void reassign_particles_() {
    for (size_t i = 0; i < N; ++i) {
      auto target_cell = ijk_to_cell_index_(pos_to_ijk_(particle_positions_0_[i]));
      if (target_cell != cell_of_particle_[i]) {
        remove_particle_from_cell_(i, cell_of_particle_[i]);
        add_particle_to_cell_(i, target_cell);
      }
    }
  }
  void get_all_in_cells_(const std::array<Vec3<size_t>, 27>& cells, std::vector<size_t>& storage) {
    for (auto cell : cells) {
      int particle = FIRST_[ijk_to_cell_index_(cell)];
      while (particle != -1) {
        storage.push_back(particle);
        particle = LIST_[particle];
      }
    }
  }
  //only call this if cell is not empty, or eat gigantic overhead
  //ijk: cell to update, neighbour_particles: vector with enough capacity to store all neighbouring particles
  //returns 2*energy from that cell
  /*
  real_t step_cell_(const Vec3<size_t> ijk, std::vector<size_t>& neighbour_particles) noexcept {
    real_t energy(0);
    int particle = FIRST_[ijk_to_cell_index_(ijk)];
    get_all_in_cells_(neighbour_cell_indices_(ijk), neighbour_particles);
    while (particle != -1) {
      energy += (particle_positions_next_[particle] - particle_positions_current_[particle]).rsq() / (4. * dt_ * dt_);
      Vec3<real_t> current_force = { 0,0,0 };
      Vec3<real_t> current_particle = particle_positions_current_[particle];
      for (auto p : neighbour_particles) {
        if (p == particle) {
          continue; //no self interaction
        }
        Vec3<real_t> distance = particle_positions_current_[p] - current_particle;
        Vec3<real_t> distance_energy = particle_positions_previous_[p] - particle_positions_previous_[particle];
        periodize_distance_(distance);
        periodize_distance_(distance_energy);
        current_force += force_(distance);
        energy += potential_(distance_energy.r());
      }
      particle_positions_next_[particle] = real_t(2) * current_particle - particle_positions_previous_[particle] + dt_ * dt_ * current_force;
      particle = LIST_[particle];
    }
    return energy;
  }*/
  //Substep 1: calculate r(t+dt), v(t+dt/2)
  //Needs xi_0_ to be up to date (done by swapping with xi_2_)
  void substep_cell_1_(const Vec3<size_t> ijk, std::vector<size_t>& neighbour_particles) noexcept{
    int particle = FIRST_[ijk_to_cell_index_(ijk)];
    get_all_in_cells_(neighbour_cell_indices_(ijk), neighbour_particles);
    while(particle != -1){
      Vec3<real_t> current_force = {0 ,0 ,0};
      Vec3<real_t> current_particle = particle_positions_0_[particle];
      for(auto p: neighbour_particles){//calculate force on particle
        if(p == particle){
          continue;
        }
        else{
          Vec3<real_t> distance = particle_positions_0_[p] - current_particle;
          periodize_distance_(distance);
          current_force += force_(distance);
        }
      }
      particle_positions_1_[particle] = current_particle + particle_velocities_0_[particle] * dt_ + (current_force - xi_0_ * particle_velocities_0_[particle]) * dt_ * dt_ * 0.5;
      particle_velocities_1_[particle] = particle_velocities_0_[particle] + (current_force - xi_0_ * particle_velocities_0_[particle]) * 0.5 * dt_;
      particle = LIST_[particle];
    }
  }
  void compute_Ekin_0_() noexcept{
    Ekin_0_ = 0;
    for(auto v: particle_velocities_0_){
      Ekin_0_ += v.rsq();
    }
    Ekin_0_ *= 0.5;
  }
  //Needs Ekin_0_ and xi_0_ to be up to date
  void compute_xi_1_() noexcept{
    xi_1_ = xi_0_ + (2. * Ekin_0_ - (3. * N + 1.) * temperature_) * (dt_ / (2. * Q_));
  }
  void compute_Ekin_1_() noexcept{
    Ekin_1_ = 0;
    for(auto v: particle_velocities_1_){
      Ekin_1_ += v.rsq();
    }
    Ekin_1_ *= 0.5;
  }
  //Needs Ekin_1_ and xi_1_ to be up to date
  void compute_xi_2_() noexcept{
    xi_2_ = xi_1_ + ( 2. * Ekin_1_ - ( 3. * N + 1. ) * temperature_ ) * ( dt_ / ( 2. * Q_ ) );
  }
  //Substep 2: calculate v(t+dt)
  //Needs xi_2_ to be up to date
  //new positions at (t+dt) must have been swapped into particle_positions_0_
  //and the cell assignment must have been updated accordingly
  void substep_cell_2_(const Vec3<size_t> ijk, std::vector<size_t>& neighbour_particles) noexcept{
    int particle = FIRST_[ijk_to_cell_index_(ijk)];
    get_all_in_cells_(neighbour_cell_indices_(ijk), neighbour_particles);
    while(particle != -1){
      Vec3<real_t> current_force = {0 ,0 ,0};
      Vec3<real_t> current_particle = particle_positions_0_[particle];
      for(auto p: neighbour_particles){//calculate force on particle
        if(p == particle){
          continue;
        }
        else{
          Vec3<real_t> distance = particle_positions_0_[p] - current_particle;
          periodize_distance_(distance);
          current_force += force_(distance);
        }
      }
      particle_velocities_2_[particle] = ( particle_velocities_1_[particle] + 0.5 * dt_ * current_force ) * (1. / ( 1. + 0.5 * dt_ * xi_2_ ) );
      particle = LIST_[particle];
    }
  }

  Lennard_Jones_Potential<real_t> potential_;
  Lennard_Jones_Force<real_t> force_;
  real_t dt_;
  real_t L_;
  size_t NM_; //linear number of cells
  std::array<int, L* L* L> FIRST_;
  std::array<int, N> LIST_;
  std::array<size_t, N> cell_of_particle_; //-1: particle not assigned to any cell
  /*
  std::array<Vec3<real_t>, N> particle_positions_next_; //t+dt or t-dt
  std::array<Vec3<real_t>, N> particle_positions_current_; //t or t+dt
  std::array<Vec3<real_t>, N> particle_positions_previous_; //t-dt or t
  std::array<Vec3<real_t>, N> particle_velocities_; //~t
  */
  std::array<Vec3<real_t>, N> particle_positions_0_; //positions at t
  std::array<Vec3<real_t>, N> particle_positions_1_; //positions at t+dt
  std::array<Vec3<real_t>, N> particle_velocities_0_; //velocities at t
  std::array<Vec3<real_t>, N> particle_velocities_1_; //velocities at t+dt/2
  std::array<Vec3<real_t>, N> particle_velocities_2_; //velocities at t+dt
  real_t xi_0_; //xi at t
  real_t xi_1_; //xi at t+dt/2
  real_t xi_2_; //xi at t+dt
  real_t Ekin_0_; //Ekin at t
  real_t Ekin_1_; //Ekin at t+dt/2
  real_t temperature_; //target temperature
  real_t Q_; //heat bath mass
};

#endif
