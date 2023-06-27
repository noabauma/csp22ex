import numpy as np
import matplotlib.pyplot as plt

out_dir = 'outputs/'

def nn_sum(x, i, j, k):
    """
    Args:
        x: Spin configuration
        i, j, k: Indices describing the position of one spin

    Returns:
        Sum of the spins in x which are nearest neighbors of (i, j, k)
    """
    result = x[(i+1)%L, j, k] + x[(i-1)%L, j, k]
    result += x[i, (j+1)%L, k] + x[i, (j-1)%L, k]
    result += x[i, j, (k+1)%L] + x[i, j, (k-1)%L]

    return int(result)


def increase_energy(x, M, E):
    """
    Args:
        x: Spin configuration
        M: Magnetization of x
        E: Energy of x

    Returns:
        Updated x, M and E
    """
    # pick one spin at random
    i = int(L*np.random.rand())
    j = int(L*np.random.rand())
    k = int(L*np.random.rand())
    x_old = x[i, j, k]

    # flip the spin iff the energy increases
    nn = nn_sum(x, i, j, k)
    dE = 2*J*x_old*nn
    if dE > 0:
        x[i, j, k] *= -1
        M -= 2*x_old
        E += dE

    return x, M, E


def move(x, M, E, Ed):
    """
    Args:
        x: Spin configuration
        M: Magnetization of x
        E: Energy of x
        Ed: Demon energy

    Returns:
        Updated x, M and E, Ed after one Monte Carlo move
    """
    # pick one spin at random
    i = int(L*np.random.rand())
    j = int(L*np.random.rand())
    k = int(L*np.random.rand())
    x_old = x[i, j, k]

    # flip the spin according to the Creutz algorithm
    nn = nn_sum(x, i, j, k)
    dE = 2*J*x_old*nn
    if Emax >= Ed - dE and Ed - dE >= 0:
        x[i, j, k] *= -1
        M -= 2*x_old
        E += dE
        Ed -= dE

    return x, M, E, Ed


L = 6 # linear size of the system
J = 1 # coupling constant
Nsample = 1000 # number of samples
Nsubsweep = 3*L**3 # number of subsweeps (3 is the safety factor)
Esys_arr = np.linspace(-3*L**3, 0, 50) # intended total system energy
                                         
Emax = 40 # maximum energy difference 
Ed = 0 # intial demon energy

Ed_arr = np.zeros(Nsample*Nsubsweep) # demon energies
T_arr = np.zeros_like(Esys_arr) # temperatures
M_arr = np.zeros_like(Esys_arr) # average magnetizations
E_arr = np.zeros_like(Esys_arr) # average energies
M_err = np.zeros_like(Esys_arr) # standard deviations of the magnetizations
E_err = np.zeros_like(Esys_arr) # standard deviations of the energies
chi_arr = np.zeros_like(Esys_arr) # magnetic susceptibilities
cv_arr = np.zeros_like(Esys_arr) # heat capacities

for i, e in enumerate(range(Esys_arr.size)):
    # start with uniform initial configuration
    x = np.ones((L, L, L))
    M = L**3
    E = -3*J*L**3 # every lattice site contributes an energy -3J

    print('Energising the system...')
    # bring the system to the specified energy
    while Esys_arr[e] > E:
      x, M, E = increase_energy(x, M, E)
    print('System reached the specified energy', Esys_arr[e], '!')

    # initialise arrays holding magnetisation and energy densities
    M_data = np.zeros(Nsample)
    E_data = np.zeros(Nsample)
    
    print('Running the Monte Carlo simulation at energy...')
    M_data[0] = np.abs(M)/L**3 # initial magnetisation density
    E_data[0] = E/L**3 # initial energy density

    for n in range(1, Nsample):
        for t in range(Nsubsweep):
            x, M, E, Ed = move(x, M, E, Ed)
            Ed_arr[n*Nsubsweep+t] = Ed # save the demon energies

        M_data[n] = np.abs(M)/L**3
        E_data[n] = E/L**3

    M_arr[e] = np.mean(M_data) # average magnetization
    E_arr[e] = np.mean(E_data) # average energy
    M_err[e] = np.std(M_data) # standard deviation of the magnetization
    E_err[e] = np.std(E_data) # standard deviation of the energy
    

    T_arr[e] = 4./(np.log(1 + 4./np.mean(Ed_arr))) 
    print('The estimated temperature of system at energy', Esys_arr[e], 'is', T_arr[e])
    
    if i % 10 == 0:
        hist, bin_edges = np.histogram(Ed_arr, bins=4*np.arange(15), density=True)
        plt.figure()
        #plt.hist((hist, bin_edges), bins=bin_edges)
        plt.plot(bin_edges, np.exp(-1/T_arr[e] * bin_edges), label=r"$T=$"+str(T_arr[e]))
        plt.bar(
            (bin_edges[1:] + bin_edges[:-1]) * .5, hist / hist.max(),
            width=(bin_edges[1] - bin_edges[0]), color="grey")
        plt.xlabel(r'$E_d$')
        plt.ylabel(r'frequency')
        plt.legend()
        plt.savefig(out_dir+'hist_T'+str(T_arr[e])+'_L'+str(L)+'.pdf')
        #plt.show()
        
     

np.savetxt(out_dir+'T'+str(L)+'.txt', T_arr, newline='\r\n')
np.savetxt(out_dir+'M'+str(L)+'.txt', M_arr, newline='\r\n')
np.savetxt(out_dir+'E'+str(L)+'.txt', E_arr, newline='\r\n')
np.savetxt(out_dir+'M_err'+str(L)+'.txt', M_err, newline='\r\n')
np.savetxt(out_dir+'E_err'+str(L)+'.txt', E_err, newline='\r\n')
