import numpy as np


def random_site():
    # pick one spin at random
    i = int(L*np.random.rand())
    j = int(L*np.random.rand())
    k = int(L*np.random.rand())
    
    return i, j, k
    

def nn_sum(x, i, j, k):
    """
    Args:
        x: Spin configuration
        i, j, k: Lattice coordinates of the central spin

    Returns:
        Sum of the spins in x which are nearest neighbors of (i, j, k)
    """
    result = x[(i+1)%L, j, k] + x[(i-1)%L, j, k]
    result += x[i, (j+1)%L, k] + x[i, (j-1)%L, k]
    result += x[i, j, (k+1)%L] + x[i, j, (k-1)%L]

    return int(result)


def increase_energy(x, M, E):
    """
    Propose a spin flip update at a randomly chosen site
    (i,j,k) and accept it only if it increases the energy.
    
    Args:
        x: Spin configuration
        M: Magnetization of x
        E: Energy of x

    Returns:
        Updated x, M and E
    """
    
    i, j, k = random_site() # pick one spin at random

    # TODO: flip the spin at site (i,j,k) if and only if the energy increases
    # TODO: update magnetisation and energy after the flip
    x_old = x[i,j,k]
    h_nn = nn_sum(x, i, j, k)
    update_E = 2*J*x_old*h_nn
    if update_E > 0:
        x[i,j,k] *= -1
        M -= 2*x_old
        E += update_E

    return x, M, E


def move(x, M, E, Ed):
    """
    Single update of the Creutz algorithm.
    
    Args:
        x: Spin configuration
        M: Magnetization of x
        E: Energy of x
        Ed: Demon energy

    Returns:
        Updated x, M and E, Ed after one Monte Carlo move
    """
    
    i, j, k = random_site()

    # flip the spin according to the Creutz algorithm
    # TODO: compute energy difference dE due to flipping spin at site (i,j,k)
    # TODO: perform the Creutz step by comparing `dE` to the demon energy `Ed`
    # TODO: update `M`, `E` and also `Ed` if the spin flip is accepted
    x_old = x[i,j,k]
    h_nn = nn_sum(x, i, j, k)
    dE = 2*J*x_old*h_nn
    if Emax >= Ed - dE >= 0:
        x[i,j,k] *= -1
        
        M  -= 2*x_old
        E  += dE
        Ed -= dE

    return x, M, E, Ed


L = 6 # linear size of the system
N = L**3 # total system size
J = 1 # coupling constant
Nsample = 5000 # number of samples
Nsubsweep = 3*N # number of subsweeps (3 is the safety factor)
Esys_arr = np.linspace(-3*N, 0, 50) # array total energies of the microcanonical ensemble
                                         
Emax = 40 # maximum energy difference 
Ed = 4 # intial demon energy

Ed_arr = np.zeros((Nsample-1)*Nsubsweep) # demon energies
T_arr = np.zeros_like(Esys_arr) # temperatures
M_arr = np.zeros_like(Esys_arr) # average magnetizations
E_arr = np.zeros_like(Esys_arr) # average energies
M_err = np.zeros_like(Esys_arr) # standard deviations of the magnetizations
E_err = np.zeros_like(Esys_arr) # standard deviations of the energies
chi_arr = np.zeros_like(Esys_arr) # magnetic susceptibilities
cv_arr = np.zeros_like(Esys_arr) # heat capacities

for e in range(Esys_arr.size):
    # start with uniform initial configuration
    x = np.ones((L, L, L))
    M = N
    E = -3*J*N # every lattice site contributes an energy -3J

    print('Energising the system...')
    # bring the system to the specified energy
    while Esys_arr[e] > E:
      x, M, E = increase_energy(x, M, E)
    print('System reached the specified energy', Esys_arr[e], '!')

    # initialise arrays holding magnetisation and energy densities
    M_data = np.zeros(Nsample)
    E_data = np.zeros(Nsample)
    
    print('Running the Monte Carlo simulation at energy...')
    M_data[0] = np.abs(M)/N # initial magnetisation density
    E_data[0] = E/N # initial energy density

    for n in range(1, Nsample):
        for t in range(Nsubsweep):
            x, M, E, Ed = move(x, M, E, Ed)
            Ed_arr[(n-1)*Nsubsweep+t] = Ed # save the demon energies

        M_data[n] = np.abs(M)/N
        E_data[n] = E/N

    M_arr[e] = np.mean(M_data) # average magnetization
    E_arr[e] = np.mean(E_data) # average energy
    M_err[e] = np.std(M_data) # standard deviation of the magnetization
    E_err[e] = np.std(E_data) # standard deviation of the energy

	
    # TODO: estimate the temperature of the system from the statistics of demon energies using data `Ed_arr`
    beta = 0.25 * np.log(1 + 4.0/np.mean(Ed_arr))
    T_arr[e] = 1.0/beta
    
    print('The estimated temperature of system at energy ', Esys_arr[e], ' is ', T_arr[e])
     
	 
out_dir = 'outputs/'
np.savetxt(out_dir+'T'+str(L)+'.txt', T_arr, newline='\r\n')
np.savetxt(out_dir+'M'+str(L)+'.txt', M_arr, newline='\r\n')
np.savetxt(out_dir+'E'+str(L)+'.txt', E_arr, newline='\r\n')
np.savetxt(out_dir+'M_err'+str(L)+'.txt', M_err, newline='\r\n')
np.savetxt(out_dir+'E_err'+str(L)+'.txt', E_err, newline='\r\n')
