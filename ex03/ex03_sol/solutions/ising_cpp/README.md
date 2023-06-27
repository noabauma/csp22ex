# Monte Carlo simulation of the Ising model

## Code Deployment

### Compilation

To compile the `C++` code you can use the standard cmake-make build. To build with code optimization in this folder, the commands are

```
cmake -DCMAKE_BUILD_TYPE=Release .
make -j
```

### Execution

The compilation will produce a first executable called `ising`. You can run it by passing an input file as an argument. The default is `input.json`.

You can specify the following in the input file:

* "dimension": the dimensionality of the hypercubic lattice, currently 1, 2 and 3 dimensions are implemented.

* "Ls": List of linear system sizes to be simulated (*e.g.*, [16, 32]).

* "algorithm": currently Metropolis-Hastings and Creutz (microcanonical) updates are implemented.

  Specific to Metropolis-Hastings algorithm:

  * "thermalization-factor" (Metropolis-Hastings): factor multiplying the system size for determining the number of steps until equilibration.
  * "Ts" (Metropolis-Hastings): List of temperatures to be simulated (*e.g.*, [1.0, 2.0]).

  Specific to Creutz algorithm:

  * "num_energies" (Creutz): number of linearly spaced energy values of the microcanonical ensemble in the range $[-dJL^d, 0]$.

  * "E_max" (Creutz): energy capacity of the demon.

  * "initial E_d" (Creutz): initial demon energy.

* "sweep-factor": factor multiplying the system size for determining the steps until uncorrelated samples.
* "measurements": total number of samples for the observables at a given temperature.

The executable will produce in the folder `outputs` a file with the magnetization and energies for each value of beta and L requested in the input file.

The results from the Metropolis-Hastings algorithm can be plotted using the script `plot.py`. It can be run by providing the simulated dimension and the set of simulated linear system sizes as command line arguments, *e.g.*, `python3 plot.py d L1 L2 L3`.

For visualising the results from the Creutz algorithm, please use the script `analyse.py`. Its usage is analogous, *e.g.*, `python3 analyse.py d L1 L2 L3`.

### Implementation details
Storing the each lattice site in a single byte integer proved to be beneficial compared to
a bitmap (simply reling on `std::vector<bool>`) due to the associated bit manipulation cost. Mind that more refined techniques where an ensemble of 8 lattices is represented in a single byte and updated in parallel might be faster, but more cumberome to implement.

