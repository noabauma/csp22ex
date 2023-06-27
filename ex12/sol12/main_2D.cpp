#include <vec2.hpp>
#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <random>

//get normalized gap vector
Vec2<double> n(const Vec2<double>& p_0, const Vec2<double>& p_1){
  Vec2<double> n_un =  p_1 - p_0;
  return n_un / n_un.r();
}
//get local gap velocity
double vn_loc(const Vec2<double>& v_0, const Vec2<double>& v_1, const Vec2<double>& n){
  return dot(v_1 - v_0, n);
}
//get new local contact force
double Rn_loc_new(const Vec2<double>& p_0, const Vec2<double>& p_1, const Vec2<double>& v_0, const Vec2<double>& v_1, double radius, double dt){
  double d_t = (p_1 - p_0).r() - 2. * radius;
  double vloc = vn_loc(v_0, v_1, n(p_0, p_1));
  return std::max(0., - 1. / (2. * dt) * (d_t / dt + vloc));
}
//get particle forces from local gap force
std::pair<Vec2<double>, Vec2<double> > R(double Rn_loc, const Vec2<double>& n){
  return {(-n) * Rn_loc, n * Rn_loc};
}
//check if gap closes
bool collision(const Vec2<double>& p_0, const Vec2<double>& p_1, const Vec2<double>& v_0, const Vec2<double>& v_1, double radius, double dt){
  return ( (p_1 - p_0).r() - 2. * radius + vn_loc(v_0, v_1, n(p_0, p_1)) * dt ) <= 0.; //corrected. we look at where we would land if we didn't correct vn_loc
}
//check if particle collides with wall and handle it (elastically)
bool collide_with_walls(const Vec2<double>& p0, const Vec2<double>& v0, Vec2<double>& p1, Vec2<double>& v1, const double radius, const double dt, const double L){
  Vec2<double> new_pos = p0 + dt * v0; //new position if we have no collisions
  //left wall
  if((new_pos.x() < radius) && (v0.x() < 0.)){
    v1 = v0;
    v1.x() *= -1.;
    p1 = p0 + dt * v1;
    return true;
  }
  //right wall
  else if((new_pos.x() > L - radius) && (v0.x() > 0.)){
    v1 = v0;
    v1.x() *= -1.;
    p1 = p0 + dt * v1;
    return true;
  }
  //bottom wall
  else if((new_pos.y() < radius) && (v0.y() < 0.)){
    v1 = v0;
    v1.y() *= -1.;
    p1 = p0 + dt * v1;
    return true;
  }
  //bottom wall
  else if((new_pos.y() > L - radius) && (v0.y() > 0.)){
    v1 = v0;
    v1.y() *= -1.;
    p1 = p0 + dt * v1;
    return true;
  }
  else{ //no collision with walls
    return false;
  }
}

int main(){
  const size_t NUM_PARTICLES = 60; //we take 60 particles
  const double L = 100.; //extent of system
  std::vector<Vec2<double> > p0 (NUM_PARTICLES); //old positions
  std::vector<Vec2<double> > p1 (NUM_PARTICLES); //new positions
  std::vector<Vec2<double> > v0 (NUM_PARTICLES); //old velocities
  std::vector<Vec2<double> > v1 (NUM_PARTICLES); //new velocities
  std::vector<std::pair<size_t, size_t> > interactions; //who interacts
  std::vector<bool> collider (NUM_PARTICLES, false); //true if particle collides in this timestep
  interactions.reserve(NUM_PARTICLES/2); //maximum number of interactions

  double radius = 1.; //radius of particles
  double dt = 0.11; //timestep


  //randomly initialize positions and velocities
  std::mt19937_64 rng(41);
  std::uniform_real_distribution<double> dis_p (0.,L);
  std::uniform_real_distribution<double> dis_v (-10.,10.);
  for(size_t i = 0; i < NUM_PARTICLES; ++i){
    p0[i] = {dis_p(rng), dis_p(rng)};
    v0[i] = {dis_v(rng), dis_v(rng)};
  }

  double max_time = 150.;
  double t = 0.; //current simulation time
  while(t <= max_time){
    t += dt;
    //step 1: check what collisions occur
    for(size_t i = 0; i < p0.size(); ++i){
      if(collider[i]){ //particle was already found to collide
        continue;
      }
      bool has_collided = false;
      //handle collisions with walls, including corresponding update.
      if(collide_with_walls(p0[i], v0[i], p1[i], v1[i], radius, dt, L)){
        continue; //particle collided with wall and was updated
      }
      for(size_t j = i+1; j < p0.size(); ++j){
        if(collision(p0[i], p0[j], v0[i], v0[j], radius, dt)){ //collision between i and j
          collider[i] = true;
          collider[j] = true;
          interactions.push_back({i, j});
          has_collided = true;
          break; //only two-particle collisions allowed
        }
      }
      if(!has_collided){ //update if we have no collision
        v1[i] = v0[i];
        p1[i] = p0[i] + v1[i] * dt;
      }
    }
    //step 2: update colliding pairs. this could also have been done above.
    for(auto ipair: interactions){
      size_t i = ipair.first;
      size_t j = ipair.second;
      auto R_force = R(Rn_loc_new(p0[i], p0[j], v0[i], v0[j], radius, dt), n(p0[i], p0[j]));
      v1[i] = v0[i] + dt * R_force.first;
      p1[i] = p0[i] + v1[i] * dt;
      v1[j] = v0[j] + dt * R_force.second;
      p1[j] = p0[j] + v1[j] * dt;
    }
    //clear interactions
    interactions.clear();
    for(auto c: collider){
      c = false;
    }
    //swap old and new
    std::swap(p0, p1);
    std::swap(v0, v1);
    //print the data
    //layout: Ekin p0.x p0.y v0 p1.x p1.y v1 ... pN.x pN.y vN
    std::cout << std::accumulate(v0.begin(), v0.end(), 0., [](const double& sum, const Vec2<double>& v){return sum + 0.5 * v.rsq();}) << " ";
    for(size_t i = 0; i < p0.size(); ++i){
      std::cout << p0[i].x() << " " << p0[i].y() << " " << v0[i].r() << " ";
    }
    std::cout << std::endl;
  }
}
