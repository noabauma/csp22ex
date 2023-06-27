CXX = /usr/local/bin/g++-11

all: main

main: main.o Ising.o
	$(CXX) main.o Ising.o -o main

main.o: main.cpp Ising.cpp Ising.hpp
	$(CXX) -c main.cpp -O3 -mtune=native -march=native -std=c++2a

Ising.o: Ising.cpp Ising.hpp
	$(CXX) -c Ising.cpp -O3 -mtune=native -march=native -std=c++2a

clean:
	rm -rf main *.o *.txt
