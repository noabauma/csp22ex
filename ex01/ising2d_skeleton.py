import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

L = 10 # linear size of the system
N = L**2
N_inv = 1./N
J = 1 # coupling constant
beta = np.linspace(0.1, 0.8, 10) # inverse temperatures
Nthermalization = int(10e5) # number of thermalization steps
Nsample = 5000 # number of samples (= size of the Markov chain)
Nsubsweep = 10*N # number of subsweeps (to generate better samples)

M_arr = np.zeros_like(beta) # average magnetizations
E_arr = np.zeros_like(beta) # average energies
M_err = np.zeros_like(beta) # standard deviations of the magnetizations
E_err = np.zeros_like(beta) # standard deviations of the energies
chi_arr = np.zeros_like(beta) # magnetic susceptibilities
cv_arr = np.zeros_like(beta) # heat capacities

def nn_sum(x, i, j):
    """
    Args:
        x: Spin configuration
        i, j: Indices describing the position of one spin

    Returns:
        Sum of the spins in x which are nearest neighbors of (i, j)
    """
    s = x[i,j]
    result = s*(x[(i+1)%L,j] + x[(i-1)%L,j] + x[i,(j+1)%L] + x[i,(j-1)%L])

    return int(result)


def total_energy(x):
    """
    Args:
        x: Spin configuration

    Returns:
        Total energy of configuration x.
    """
    energy = 0
    for i in range(L):
        for j in range(L):
            # TODO: compute energy of site (i,j)
            energy += nn_sum(x,i,j)
    energy *= -J
			
    return energy


def move(x, M, E, beta_t):
    """
    Args:
        x: Spin configuration
        M: Magnetization of x
        E: Energy of x

    Returns:
        Updated x, M and E after one Monte Carlo move
    """
    # Probability look-up tables
	# TODO: optionally use probability lookup tables
	
    # TODO: pick one site at random
    i,j = np.random.randint(L,size=2)

    # Flip the spin of that site according to the Metropolis algorithm
    # TODO: compute the local magnetic field at site (i,j) due to nearest-neighbours
    dE = 2*J*nn_sum(x,i,j)

	# TODO: Compute the Metropolis acceptance probability `R` and compare it to a random number in [0,1)
    R = np.exp(-beta_t*dE)

    eta = np.random.rand()
    if R > eta:
        # TODO: flip the spin
        x[i,j] *= -1

        # TODO: update the magnetisation and energy
        #inefficient way
        M = N_inv*np.sum(x)
        E = total_energy(x)

        #cool way
        #M += 2*N_inv*x[i,j] the cool way
        #E = total_energy(x)

    return x, M, E



def main():
    # calculate the relevant physical quantities for different temperatures
    for t in range(beta.size):
        beta_t = beta[t]
        print('Running for inverse temperature =', beta_t)

        # TODO: Generate a random initial configuration
        # TODO: compute its magnetisation and energy
        x = np.random.randint(2, size=(L,L))
        x[x==0] = -1

        # TODO: run a thermalisation loop
        M = N_inv*np.sum(x)
        E = total_energy(x)

        for Nt in tqdm(range(Nthermalization)):
            x, M, E = move(x, M, E, beta_t)

        # measurement of M and E

        print('Sampling M and E ...')

        M_data = np.zeros(Nsample)
        E_data = np.zeros(Nsample)

        M_data[0] = np.abs(M)/N
        E_data[0] = E/N

        for n in tqdm(range(1, Nsample)):
            for N_ in range(Nsubsweep):     #intended "N"?
                x, M, E = move(x, M, E, beta_t)

            M_data[n] = np.abs(M)/N
            E_data[n] = E/N

        M_arr[t] = np.mean(M_data) # average magnetization
        E_arr[t] = np.mean(E_data) # average energy
        # TODO: use the fluctuation dissipation to compute the 
        # specific heat and susceptibility from the M and E data
        M_err[t] = np.std(M_data)
        E_err[t] = np.std(E_data)
        
        chi_arr[t] = beta_t*beta_t/(N*(np.mean(M_data**2)-M_err[t]**2))
        cv_arr[t] = beta_t*beta_t/(N*(np.mean(E_data**2)-E_err[t]**2))
        #cv_arr[t] = beta_t*beta_t/(N*E_err[t])
        

    # plot magnetization
    plt.figure()
    plt.errorbar(beta, M_arr, yerr=M_err)
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$M$')
    plt.savefig('M.png')

    # TODO: plot magnetic susceptibility
    plt.figure()
    plt.plot(beta, chi_arr)
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\chi$')
    plt.savefig('chi.png')

    # TODO: plot energy
    plt.figure()
    plt.errorbar(beta, E_arr, yerr=E_err)
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$E$')
    plt.savefig('E.png')

    # TODO: plot heat capacity
    plt.figure()
    plt.plot(beta, cv_arr)
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$C_v$')
    plt.savefig('cv.png')


if __name__ == "__main__":
    main()