CXX=/usr/local/bin/g++-11
all: sanity_check speed bonus

sanity_check: sanity_check.o Ising.o
	$(CXX) sanity_check.o Ising.o -o sanity_check -Wall -Wpedantic -O3 -mtune=native -march=native -ffast-math -fstrict-aliasing -funroll-loops

sanity_check.o: sanity_check.cpp
	$(CXX) -c sanity_check.cpp -std=c++20 -Wall -Wpedantic -O3 -mtune=native -march=native -ffast-math -fstrict-aliasing -funroll-loops

speed: speed.o Ising.o
		$(CXX) speed.o Ising.o -o speed -Wall -Wpedantic -O3 -mtune=native -march=native -ffast-math -fstrict-aliasing -funroll-loops

speed.o: speed.cpp
		$(CXX) -c speed.cpp -std=c++20 -Wall -Wpedantic -O3 -mtune=native -march=native -ffast-math -fstrict-aliasing -funroll-loops

bonus: bonus.o Ising.o
	$(CXX) bonus.o Ising.o -o bonus -Wall -Wpedantic -O3 -mtune=native -march=native -ffast-math -fstrict-aliasing -funroll-loops

bonus.o: bonus.cpp
		$(CXX) -c bonus.cpp -std=c++20 -Wall -Wpedantic -O3 -mtune=native -march=native -ffast-math -fstrict-aliasing -funroll-loops

Ising.o: Ising.cpp Ising.hpp
	$(CXX) -c Ising.cpp -std=c++20 -Wall -Wpedantic -O3 -mtune=native -march=native -ffast-math -fstrict-aliasing -funroll-loops

clean:
	rm -rf sanity_check speed bonus *.o

cleanall:
	rm -rf sanity_check speed bonus *.o *.csv
