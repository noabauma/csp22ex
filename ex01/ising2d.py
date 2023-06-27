import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import numba
import platform

L = 10 # linear size of the system
N = L**2
N_inv = 1./N
J = 1 # coupling constant
n_temps = 10
temps = np.linspace(0.1, 8, n_temps)
beta = 1.0/temps
#beta = np.linspace(0.1, 0.8, n_temps) # inverse temperatures
Nthermalization = int(10e5) # number of thermalization steps
Nsample = 5000 # number of samples (= size of the Markov chain)
Nsubsweep = 10*N # number of subsweeps (to generate better samples)

M_arr = np.zeros_like(beta) # average magnetizations
E_arr = np.zeros_like(beta) # average energies
M_err = np.zeros_like(beta) # standard deviations of the magnetizations
E_err = np.zeros_like(beta) # standard deviations of the energies
chi_arr = np.zeros_like(beta) # magnetic susceptibilities
cv_arr = np.zeros_like(beta) # heat capacities
chi_arr2 = np.zeros_like(beta) # magnetic susceptibilities
cv_arr2 = np.zeros_like(beta) # heat capacities

#calculate the look-up table for the acceptance probability R
R_mat = np.empty((n_temps,5), dtype=float)
for i in range(n_temps):
    for j in range(5):
        R_mat[i,j] = np.exp(-beta[i]*2*J*(2*j-4))

@numba.njit()
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


@numba.njit()
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


@numba.jit()
def move(x, M, E, beta_t, t):
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
    nn = nn_sum(x,i,j)
    dE = 2*J*nn

    if dE <= 0.0:
        # TODO: flip the spin
        E += J*(nn_sum(x,i,j) + nn_sum(x,(i+1)%L,j) + nn_sum(x,(i-1)%L,j) + nn_sum(x,i,(j+1)%L) + nn_sum(x,i,(j-1)%L))
        x[i,j] *= -1

        # TODO: update the magnetisation and energy

        #cool way
        M += 2*N_inv*x[i,j]
        E += -J*(nn_sum(x,i,j) + nn_sum(x,(i+1)%L,j) + nn_sum(x,(i-1)%L,j) + nn_sum(x,i,(j+1)%L) + nn_sum(x,i,(j-1)%L))

    else:
        # TODO: Compute the Metropolis acceptance probability `R` and compare it to a random number in [0,1)
        #R = np.exp(-beta_t*dE)
        R = R_mat[t,int(0.5*(nn+4))]

        eta = np.random.rand()
        if R > eta:
            # TODO: flip the spin
            E += J*(nn_sum(x,i,j) + nn_sum(x,(i+1)%L,j) + nn_sum(x,(i-1)%L,j) + nn_sum(x,i,(j+1)%L) + nn_sum(x,i,(j-1)%L))
            x[i,j] *= -1

            # TODO: update the magnetisation and energy

            #cool way
            M += 2*N_inv*x[i,j]
            E += -J*(nn_sum(x,i,j) + nn_sum(x,(i+1)%L,j) + nn_sum(x,(i-1)%L,j) + nn_sum(x,i,(j+1)%L) + nn_sum(x,i,(j-1)%L))

    return x, M, E


@numba.jit(parallel=True)
def main():
    # calculate the relevant physical quantities for different temperatures
    for t in numba.prange(beta.size):
        beta_t = beta[t]
        print('Running for inverse temperature =', beta_t)

        # TODO: Generate a random initial configuration
        # TODO: compute its magnetisation and energy
        x = np.random.randint(2, size=(L,L))
        x[x==0] = -1

        M = np.sum(x)
        E = total_energy(x)

        # TODO: run a thermalisation loop

        for Nt in tqdm(range(Nthermalization)):
            x, M, E = move(x, M, E, beta_t, t)

        # measurement of M and E

        print('Sampling M and E ...')

        M_data = np.zeros(Nsample)
        E_data = np.zeros(Nsample)

        M_data[0] = np.abs(M)/N
        E_data[0] = E/N

        for n in tqdm(range(1, Nsample)):
        #for n in range(1, Nsample):
            for N_ in range(Nsubsweep):     #intended "N"?
                x, M, E = move(x, M, E, beta_t, t)

            M_data[n] = np.abs(M)/N
            E_data[n] = E/N

        M_arr[t] = np.mean(M_data) # average magnetization
        E_arr[t] = np.mean(E_data) # average energy
        # TODO: use the fluctuation dissipation to compute the 
        # specific heat and susceptibility from the M and E data
        M_err[t] = np.std(M_data)
        E_err[t] = np.std(E_data)
        
        chi_arr[t] = N*beta_t*(np.mean(M_data**2)-M_arr[t]**2)
        chi_arr2[t] = N*beta_t*M_err[t]
        cv_arr[t] = beta_t*beta_t*(np.mean(E_data**2)-E_arr[t]**2)
        cv_arr2[t] = beta_t*beta_t*E_err[t]
        
    #temps = 1.0/beta

    fig, axs = plt.subplots(2,3)
    fig.suptitle('plots')

    axs[0,0].errorbar(temps, M_arr, yerr=M_err)
    axs[0,0].set_xlabel(r'$T$')
    axs[0,0].set_ylabel(r'$\langle |M| \rangle$')
    axs[0,0].set_title('magnetization')

    axs[0,1].plot(temps, chi_arr)
    axs[0,1].set_xlabel(r'$T$')
    axs[0,1].set_ylabel(r'$\chi$')
    axs[0,1].set_title('magnetic susceptibility')

    axs[0,2].plot(temps, chi_arr2)
    axs[0,2].set_xlabel(r'$T$')
    axs[0,2].set_ylabel(r'$\chi$')
    axs[0,2].set_title('magnetic susceptibility 2')

    axs[1,0].errorbar(temps, E_arr, yerr=E_err)
    axs[1,0].set_xlabel(r'$T$')
    axs[1,0].set_ylabel(r'$\langle E \rangle$')
    axs[1,0].set_title('energy')

    axs[1,1].plot(temps, cv_arr)
    axs[1,1].set_xlabel(r'$T$')
    axs[1,1].set_ylabel(r'$C_v$')
    axs[1,1].set_title('heat capacity')

    axs[1,2].plot(temps, cv_arr2)
    axs[1,2].set_xlabel(r'$T$')
    axs[1,2].set_ylabel(r'$C_v$')
    axs[1,2].set_title('heat capacity 2')

    if "WSL" in platform.uname().release:
        plt.savefig("plots.png")
    else:
        plt.show()


if __name__ == "__main__":
    main()