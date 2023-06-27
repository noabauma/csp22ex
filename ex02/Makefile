CXX = g++
CPPFLAGS = -O3 -Isrc/ -I/usr/include/eigen3 -march=native -ffast-math -funroll-loops \
		 -Wall -Wextra \

.PHONY: clean

all: task1 task2_3 task4

task1: ising3d_task1.cpp
	$(CXX) $(CPPFLAGS) -o ising3d_task1 ising3d_task1.cpp

task2_3: ising3d_task2_3.cpp
	$(CXX) $(CPPFLAGS) -o ising3d_task2_3 ising3d_task2_3.cpp

task4: ising3d_task4.cpp
	$(CXX) $(CPPFLAGS) -o ising3d_task4 ising3d_task4.cpp 

run_task1: ising3d_task1
	./ising3d_task1 4.1
	python3 plot_task1.py task1.csv 4.1

run_task2_3: ising3d_task2_3
	./ising3d_task2_3 4.1
	
run_task4: ising3d_task4
	./ising3d_task4 4.1 8

clean:
	rm -f ising3d_task1 ising3d_task2_3 ising3d_task4
