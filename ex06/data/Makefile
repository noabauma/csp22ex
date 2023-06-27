CXX = g++
CPPFLAGS = -O3 -I/usr/include/eigen3 \
		 -Wall -Wextra \

.PHONY: clean

main: ising2d.cpp
	$(CXX) $(CPPFLAGS) -o ising2d ising2d.cpp 

run: ising2d
	./ising2d 2.1 1
	
plot:
	python3 plot.py data.csv
	
clean:
	rm -f ising2d
