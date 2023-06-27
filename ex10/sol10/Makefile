CXX = /usr/local/bin/g++-11

INCLUDE = -Iinclude
FLAGS = -O3 -mtune=native -march=native -fstrict-aliasing -finline-functions

all: main

main: main.cpp include/lennard_jones.hpp include/md_box_nh.hpp include/vec3.hpp
	${CXX} main.cpp -o main ${INCLUDE} ${FLAGS}

clean:
	rm -rf *.o main *.txt
