cmake_minimum_required(VERSION 3.0)
project(heisenberg3D)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native -ftree-vectorize")

include_directories(include)
add_library(src src/heisenberg_lattice.cpp src/json11.cpp src/metropolis.cpp src/wolff.cpp
            src/time_series.cpp src/clock.cpp)

add_executable(heisenberg3d main_heisenberg.cpp)
target_link_libraries(heisenberg3d src)

add_executable(autocorrelation main_autocorrelation.cpp)
target_link_libraries(autocorrelation src)
