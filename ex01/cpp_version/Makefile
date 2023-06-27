CXX = g++
CPPFLAGS = -O3 -I/usr/include/eigen3 -march=native -ffast-math -funroll-loops \
		 -Wall -Wextra \

.PHONY: clean

main: ising2d.cpp
	$(CXX) $(CPPFLAGS) -o main ising2d.cpp 

main2: ising2d_v2.cpp
	$(CXX) $(CPPFLAGS) -o main2 ising2d_v2.cpp 

run: main
	./main 2.1

run2: main2
	./main2 2.1
	
plot:
	python3 plot.py data.csv
	
clean:
	rm -f main main2
