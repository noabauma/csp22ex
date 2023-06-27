CXX=/usr/local/bin/g++-11
all: main correlation_analysis

correlation_analysis: correlation_analysis.o Ising.o
	$(CXX) correlation_analysis.o Ising.o -o correlation_analysis -Wall -Wpedantic -O3 -mtune=native -march=native

main: main.o Ising.o
	$(CXX) main.o Ising.o -o main -Wall -Wpedantic -O3 -mtune=native -march=native

main.o: main.cpp
	$(CXX) -c main.cpp -std=c++20 -Wall -Wpedantic -O3 -mtune=native -march=native

Ising.o: Ising.cpp Ising.hpp
	$(CXX) -c Ising.cpp -std=c++20 -Wall -Wpedantic -O3 -mtune=native -march=native

clean:
	rm -rf main correlation_analysis *.o *.csv
