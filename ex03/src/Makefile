CXX = g++
CPPFLAGS = -O3 -march=native -ffast-math -funroll-loops \
		 -Wall -Wextra \

.PHONY: clean

main: lattice.cpp
	$(CXX) $(CPPFLAGS) -o lattice lattice.cpp 

test: tensor_test.cpp
	$(CXX) $(CPPFLAGS) -I/usr/include/eigen-master -o test tensor_test.cpp 


run: lattice
	./lattice 
	
clean:
	rm -f lattice
