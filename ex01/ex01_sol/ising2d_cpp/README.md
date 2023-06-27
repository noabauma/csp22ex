# Metropolis-Hastings algorithm for the 2-d Ising model on square lattice

## Code Deployment

### Compilation

To compile the `C++` code you can use the standard cmake-make build. To build with code optimization in this folder the commands are

```
cmake -DCMAKE_BUILD_TYPE=Release .
make -j
```

### Execution

The compilation will produce a first executable called `ising`. You can run it by passing an input file as an argument. The default is `input.json`.
The executable will produce in the folder `outputs` a file with the magnetization and energies for each value of beta and L requested in the input file.

The results can be plotted using the script `plot.py`. It can be run by providing the set of simulated linear system sizes as command line arguments, *e.g.*, `python3 plot.py L1 L2 L3`.

### Implementation details
Storing the each lattice site in a single byte integer proved to be beneficial compared to
a bitmap (simply reling on `std::vector<bool>`) due to the associated bit manipulation cost. Mind that more refined techniques where an ensemble of 8 lattices is represented in a single byte and updated in parallel might be faster, but more cumberome to implement.

