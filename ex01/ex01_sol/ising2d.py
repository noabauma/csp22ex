import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def nn_sum(x, i, j):
    """
    Args:
        x: Spin configuration
        i, j: Indices describing the position of one spin

    Returns:
        Sum of the spins in x which are nearest neighbors of (i, j)
    """
    result = x[(i+1)%L, j] + x[(i-1)%L, j]
    result += x[i, (j+1)%L] + x[i, (j-1)%L]
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
            energy += -J*x[i,j]*nn_sum(x, i, j)/2
    return energy


def move(x, M, E):
    """
    Args:
        x: Spin configuration
        M: Magnetization of x
        E: Energy of x

    Returns:
        Updated x, M and E after one Monte Carlo move
    """
    # pick one spin at random
    i = int(L*np.random.rand())
    j = int(L*np.random.rand())
    x_old = x[i, j]

    # flip the spin according to the Metropolis algorithm
    h_nn = nn_sum(x, i, j) # local magnetic field due to nearest-neighbours
    if x_old == 1:
        R = table_spin_up[int((h_nn+4)/2)] # Metropolis acceptance probability
    else:
        R = table_spin_down[int((h_nn+4)/2)]
    eta = np.random.rand()
    if R > eta:
        x[i, j] *= -1
        M -= 2*x_old
        E += 2*J*x_old*h_nn

    return x, M, E


L = 10 # linear size of the system
N = L**2
J = 1 # coupling constant
betas = np.linspace(0.1, 0.8, 10) # inverse temperatures
#betas = [1.0/2.1]
Nthermalization = int(10e5) # number of thermalization steps
Nsample = 5000 # number of samples (= size of the Markov chain)
Nsubsweep = 10*N # number of subsweeps (to generate better samples)

M_arr = np.zeros_like(betas) # average magnetizations
E_arr = np.zeros_like(betas) # average energies
M_err = np.zeros_like(betas) # standard deviations of the magnetizations
E_err = np.zeros_like(betas) # standard deviations of the energies
chi_arr = np.zeros_like(betas) # magnetic susceptibilities
cv_arr = np.zeros_like(betas) # heat capacities

# calculate the relevant physical quantities for different temperatures
for t, beta in enumerate(betas):
    print('Running for inverse temperature =', beta)

    # Probability look-up tables
    table_spin_up = np.exp(-2*J*beta*np.array([-4, -2, 0, 2, 4]))
    table_spin_down = np.exp(+2*J*beta*np.array([-4, -2, 0, 2, 4]))

    # random initial configuration
    #x = np.random.randint(2, size=(L, L))*2 - 1 
    #M = np.sum(x)
    #E = total_energy(x)

    x = np.ones((L, L)) # initialise with uniform configuration
    M = N # magnetisation for uniform configuration
    E = -2*J*N # every lattice site contributes an energy -2J
    for i in range(L):
        for j in range(L):
            if 0.5 > np.random.rand():
                x[i, j] = -1
                M -= 2 
                E += 2*J*nn_sum(x, i, j)

    print('Thermalizing ...')

    # thermalization
    for nt in tqdm(range(Nthermalization)):
        x, M, E = move(x, M, E)

    print('Sampling M and E ...')

    # measurement of M and E
    M_data = np.zeros(Nsample)
    E_data = np.zeros(Nsample)

    M_data[0] = np.abs(M)/N
    E_data[0] = E/N

    for n in tqdm(range(1, Nsample)):
        for _ in range(Nsubsweep):
            x, M, E = move(x, M, E)
        M_data[n] = np.abs(M)/N
        E_data[n] = E/N

    M_arr[t] = np.mean(M_data) # average magnetization
    E_arr[t] = np.mean(E_data) # average energy
    M_err[t] = np.std(M_data) # standard deviation of the magnetization
    E_err[t] = np.std(E_data) # standard deviation of the energy
    chi_arr[t] = N*beta*M_err[t]**2 # magnetic susceptibility per spin
    cv_arr[t] = beta**2*E_err[t]**2 # heat capacity per spin
    
    #print(M_arr, M_err, chi_arr, E_arr, E_err, cv_arr)

# plot magnetization
plt.figure()
plt.errorbar(betas, M_arr, yerr=M_err)
plt.xlabel('$\\beta$')
plt.ylabel('$M$')
plt.savefig('M.pdf')

# plot magnetic susceptibility
plt.figure()
plt.plot(betas, chi_arr, '.-')
plt.xlabel('$\\beta$')
plt.ylabel('$\\chi$')
plt.savefig('chi.pdf')

# plot energy
plt.figure()
plt.errorbar(betas, E_arr, yerr=E_err)
plt.xlabel('$\\beta$')
plt.ylabel('$E$')
plt.savefig('E.pdf')

# plot heat capacity
plt.figure()
plt.plot(betas, cv_arr, '.-')
plt.xlabel('$\\beta$')
plt.ylabel('$c_V$')
plt.savefig('cv.pdf')
