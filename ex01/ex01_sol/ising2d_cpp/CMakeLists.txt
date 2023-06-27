cmake_minimum_required(VERSION 3.0)
project(ising2d)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native -ftree-vectorize")

include_directories(include)
add_library(src src/ising_lattice.cpp src/metropolis.cpp src/json11.cpp)

add_executable(ising main_ising.cpp)
target_link_libraries(ising src)
