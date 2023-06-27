CXX = g++
CPPFLAGS = -O3 -Isrc/ -I/usr/include/eigen3 -march=native -ffast-math -funroll-loops \
		 -Wall -Wextra \

.PHONY: clean

all: ising3D

main: ising3D_creutz.cpp
	$(CXX) $(CPPFLAGS) -o ising3D_creutz ising3D_creutz.cpp

run_task1: ising3D_creutz
	./ising3D_creutz -300.0

run_task2: ising3D_creutz
	source run_task2.sh

clean:
	rm -f ising3D
