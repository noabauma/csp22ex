CXX = g++
CPPFLAGS = -O3 -march=native -ffast-math -funroll-loops \
		 -Wall -Wextra \

.PHONY: clean

main: lattice.cpp
	$(CXX) $(CPPFLAGS) -o lattice lattice.cpp 


run: lattice
	./lattice 
	
clean:
	rm -f lattice
