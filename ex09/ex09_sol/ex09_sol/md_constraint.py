import os
import vtktools
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from tqdm import tqdm

data_dir = "simu"
if not os.path.exists(data_dir):
    os.makedirs(data_dir)


def potential(ri, rj):
    """
    Args:
        ri: Position vector of particle i
        rj: Position vector of particle j (i must not be equal to j)

    Returns:
        (Lennard-Jones) potential energy between both particles
    """
    r_vec = ri - rj

    # account for periodic boundary conditions
    for k in range(3):
        r_k = ri[k]-rj[k]
        if abs(r_k) > L/2.0:
            r_vec[k] = -np.sign(r_k)*(L - abs(r_k))

    r = np.linalg.norm(r_vec)

    return 4*((1/r)**12 - (1/r)**6)


def energy(r_current, v_current):
    """
    Args:
        r_current: Current particle positions
        v_current: Current particle velocities

    Returns:
        Total energy of the system
    """
    E_kin = 0.5*sum(sum(np.square(v_current)))

    E_pot = 0
    V_c = 4*((1/rc)**12-(1/rc)**6) # cutoff potential

    for i in range(len(r_current)):
        for j in range(len(r_current) - i - 1):
            E_pot += potential(r_current[i, :], r_current[i+j+1, :]) - V_c

    return E_kin + E_pot
    
    
def pot_energy(r_current1, r_current2):
    """
    Args:
        r_current: Current particle positions

    Returns:
        Total potential energy of the system
    """
    
    E_pot = 0
    V_c = 4*((1/rc)**12-(1/rc)**6) # cutoff potential

    for i in range(N):
        for j in range(N):
            if i != j:
                E_pot += potential(r_current1[i, :], r_current1[j, :]) - V_c
                E_pot += potential(r_current1[i, :], r_current2[j, :]) - V_c
                
                E_pot += potential(r_current2[i, :], r_current1[j, :]) - V_c
                E_pot += potential(r_current2[i, :], r_current2[j, :]) - V_c

    return E_pot


def multiplier(r, r_tilde):
    """
    Args:
        r: Particle positions before Verlet step part 1
        r_tilde: Particle positions after Verlet step part 1

    Returns:
        Lagrange multiplier
    """
    dot1 = np.dot(r, r_tilde)
    dot2 = np.dot(r, r)
    dot3 = np.dot(r_tilde, r_tilde)

    term1 = -dot1
    term2 = np.sqrt(dot1**2 - d**2*(dot3 - d**2))

    if np.abs(term1 + term2) <= np.abs(term1 - term2):
        result = 1/(2*d**2)*(term1 + term2)
    else:
        result = 1/(2*d**2)*(term1 - term2)

    return result


def stepVerlet(r_previous1, r_previous2, r_current1, r_current2):
    """
    Args:
        r_previous1: Particle positions at time t-dt (atom 1)
        r_previous2: Particle positions at time t-dt (atom 2)
        r_current1: Particle positions at time t (atom 1)
        r_current2: Particle positions at time t (atom 2)

    Returns:
        Updated positions as well as velocities and forces according to the
        Verlet scheme
    """
    # compute the forces acting on every particle (=acceleration)
    r_current1_pbc = r_current1%L
    r_current2_pbc = r_current2%L

    F1 = np.zeros((N, 3))
    F2 = np.zeros((N, 3))

    for i in range(N):
        for j in range(N):
            if i != j:
                F1[i, :] += force(r_current1_pbc[i, :], r_current1_pbc[j, :])
                F1[i, :] += force(r_current1_pbc[i, :], r_current2_pbc[j, :])
                F2[i, :] += force(r_current2_pbc[i, :], r_current2_pbc[j, :])
                F2[i, :] += force(r_current2_pbc[i, :], r_current1_pbc[j, :])

    r_next1 = np.zeros((N, 3))
    r_next2 = np.zeros((N, 3))
    del_r1 = np.zeros((N, 3))
    del_r2 = np.zeros((N, 3))

    # compute the new positions using the Verlet scheme (intermolecular part)
    for i in range(N):
        r_next1[i, :] = 2*r_current1[i, :] - r_previous1[i, :] + F1[i, :]*dt**2
        r_next2[i, :] = 2*r_current2[i, :] - r_previous2[i, :] + F2[i, :]*dt**2

    # compute the new positions using the Verlet scheme (intramolecular part)
    for i in range(N):
        r = r_current1_pbc[i, :] - r_current2_pbc[i, :]
        r_tilde = r_next1[i, :]%L - r_next2[i, :]%L
        # account for periodic boundary conditions
        for j in range(3):
            if np.abs(r[j]) > L/2.0:
                r[j] = -np.sign(r[j])*(L - np.abs(r[j]))
            if np.abs(r_tilde[j]) > L/2.0:
                r_tilde[j] = -np.sign(r_tilde[j])*(L - np.abs(r_tilde[j]))

        lam = multiplier(r, r_tilde)
        r_next1[i, :] = r_next1[i, :] + lam*r
        r_next2[i, :] = r_next2[i, :] - lam*r

        del_r1[i, :] = r_next1[i, :] - r_previous1[i, :]
        del_r2[i, :] = r_next2[i, :] - r_previous2[i, :]

        if any(r_current1[i,:] != r_current1_pbc[i,:]):
            # particle i went across the boundary in the previous time step
            r_current1[i,:] = r_current1_pbc[i,:]
            r_next1[i,:] = r_next1[i,:]%L
        if any(r_current2[i,:] != r_current2_pbc[i,:]):
            r_current2[i,:] = r_current2_pbc[i,:]
            r_next2[i,:] = r_next2[i,:]%L

    v_current1 = del_r1/(2*dt)
    v_current2 = del_r2/(2*dt)

    return r_current1, r_current2, v_current1, v_current2, r_next1, r_next2, F1, F2


def force(ri, rj):
    """
    Args:
        ri: Position vector of particle i
        rj: Position vector of particle j (i must not be equal to j)

    Returns:
        (Lennard-Jones) force between both particles
    """
    r_vec = ri - rj

    # account for periodic boundary conditions
    for k in range(3):
        r_k = ri[k]-rj[k]
        if abs(r_k) > L/2.0:
            r_vec[k] = -np.sign(r_k)*(L - abs(r_k))

    r = np.linalg.norm(r_vec)

    if r > rc:
        return np.zeros(3)
    else:
        return 24*(2*(1/r)**14 - (1/r)**8)*r_vec


"""
Parameters
"""
N = 60 # number of molecules
d = 0.2 # molecular distance
#L = 20 # box length
rc = 2.5 # cutoff-length
l = 4 # number of unit cells in one dimension
L_unit = 1 * rc # length unit
L = l * L_unit # linear box length

T = 1000 # simulation steps
dt = 10**-3 # time step

"""
Initialization
"""
potenergy_arr = np.zeros(T)
kinenergy_arr = np.zeros(T)

# initialise the particles carefully so that they do not overlap
all_cell_indices = np.array([[i, j, k] for i in range(l) for j in range(l) for k in range(l)])
cell_indices = all_cell_indices[np.random.choice(all_cell_indices.shape[0], N, replace=False), :]

r_current1 = np.zeros((N, 3)) # particle positions at time t0
for i in range(N):
    #r_current1[i] = np.random.rand(3)*L
    r_current1[i] = np.random.rand(3) + cell_indices[i] * rc

r_current2 = deepcopy(r_current1)
r_current2[:, 0] = (r_current2[:, 0] - d) % L # second atom of the molecule is close to the first atom

v_current = np.random.normal(0, 5, (N, 3)) # particle velocities at time t0
# particle positions at time t0+dt
r_next1 = r_current1 + v_current*dt
r_next2 = r_current2 + v_current*dt

"""
Simulation
"""
vtk_writer1 = vtktools.VTK_XML_Serial_Unstructured()
vtk_writer2 = vtktools.VTK_XML_Serial_Unstructured()
for t in tqdm(range(T)):
    r_current1, r_current2, v_current1, v_current2, r_next1, r_next2, \
        F_ij1, F_ij2 = stepVerlet(r_current1, r_current2, r_next1, r_next2)

    potenergy_arr[t] = 0.5 * pot_energy(r_current1%L, r_current2%L)
    kinenergy_arr[t] = 0.5 *sum(sum(np.square(v_current1))) + 0.5 * sum(sum(np.square(v_current2)))

    r_x1 = r_current1[:, 0]
    r_y1 = r_current1[:, 1]
    r_z1 = r_current1[:, 2]
    F_x1 = F_ij1[:, 0]
    F_y1 = F_ij1[:, 1]
    F_z1 = F_ij1[:, 2]
    r_x2 = r_current2[:, 0]
    r_y2 = r_current2[:, 1]
    r_z2 = r_current2[:, 2]
    F_x2 = F_ij2[:, 0]
    F_y2 = F_ij2[:, 1]
    F_z2 = F_ij2[:, 2]
    vtk_writer1.snapshot("simu/MD"+str(t)+".vtu", r_x1, r_y1, r_z1, \
        x_force=F_x1, y_force=F_y1, z_force=F_z1)
    vtk_writer2.snapshot("simu/MD"+str(t)+"n.vtu", r_x2, r_y2, r_z2, \
        x_force=F_x2, y_force=F_y2, z_force=F_z2)

vtk_writer1.writePVD("MD.pvd")
vtk_writer2.writePVD("MD1.pvd")

"""
Plotting the system energy
"""
energy_arr = potenergy_arr + kinenergy_arr

plt.figure()
#plt.semilogy(energy_arr)
plt.plot(energy_arr)
plt.ylim(0, 1.1*np.max(energy_arr))
plt.xlabel('Timesteps')
plt.ylabel('Energy')

plt.figure()
plt.plot(potenergy_arr, color='green')
plt.xlabel('Timesteps')
plt.ylabel('Potential Energy')

plt.figure()
plt.plot(kinenergy_arr, color='red')
plt.xlabel('Timesteps')
plt.ylabel('Kinetic Energy')

plt.show()
