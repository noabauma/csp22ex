CXX = g++
CPPFLAGS = -O3 -Isrc/ -I/usr/include/eigen3 -march=native -ffast-math -funroll-loops \
		 -Wall -Wextra \

.PHONY: clean

wolff: ising3d_wolff.cpp
	$(CXX) $(CPPFLAGS) -o ising3d_wolff ising3d_wolff.cpp

mr2t2: ising3d_mr2t2.cpp
	$(CXX) $(CPPFLAGS) -o ising3d_mr2t2 ising3d_mr2t2.cpp

wolff_task3: ising3d_wolff_task3.cpp
	$(CXX) $(CPPFLAGS) -o ising3d_wolff_task3 ising3d_wolff_task3.cpp

mr2t2_task3: ising3d_mr2t2_task3.cpp
	$(CXX) $(CPPFLAGS) -o ising3d_mr2t2_task3 ising3d_mr2t2_task3.cpp

run_task2: ising3d_wolff
	source run_task2.sh

clean:
	rm -f ising3d_wolff
