import os
import vtktools
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from tqdm import tqdm


data_dir = "simu"
if not os.path.exists(data_dir):
    os.makedirs(data_dir)


def r_rel_pbc(ri, rj):
    """
    Compute the relative position of particles i and j in a box
    with periodic boundary conditions.
    
    Args:
        ri: Position vector of particle i
        rj: Position vector of particle j (i must not be equal to j)
    
    Returns: 
        the shortest relative distance with correct orientation
        between particles i and j in periodic boundary conditions (PBC)
    """
    r_vec = ri - rj # relative distance without PBC
    
    for k in range(3):
        r_k = ri[k]-rj[k]
        if abs(r_k) > L/2.0:
            r_vec[k] = -np.sign(r_k)*(L - abs(r_k))
    
    return r_vec
    

def potential(ri, rj):
    """
    Args:
        ri: Position vector of particle i
        rj: Position vector of particle j (i must not be equal to j)

    Returns:
        (Lennard-Jones) potential energy between particles
    """
    r_vec = r_rel_pbc(ri, rj)
    r = np.linalg.norm(r_vec)
    
    if r > rc: # cut-off
        return 4*((1/rc)**12 - (1/rc)**6)
    else:
        return 4*((1/r)**12 - (1/r)**6)
    
    
def force(ri, rj):
    """
    Args:
        ri: Position vector of particle i
        rj: Position vector of particle j (i must not be equal to j)

    Returns:
        (Lennard-Jones) force vector
    """
    r_vec = r_rel_pbc(ri, rj)
    r = np.linalg.norm(r_vec)

    if r > rc: # cut-off
        return np.zeros(3)
    else:
        return 24*(2*(1/r)**14 - (1/r)**8)*r_vec


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


def stepVerlet(r_previous, r_current):
    """
    Args:
        r_previous: Particle positions at time t-dt
        r_current: Particle positions at time t

    Returns:
        Updated positions as well as velocities and forces according to the
        Verlet scheme
    """
    # compute the forces acting on every particle (=acceleration)

    r_current_pbc = r_current%L

    F = np.zeros((N, 3))

    for i in range(N):
        for j in range(N):
            if i != j:
                F[i, :] += force(r_current_pbc[i, :], r_current_pbc[j, :])

    r_next = np.zeros((N, 3))
    del_r = np.zeros((N, 3))

    # compute the new positions using the Verlet scheme
    for i in range(N):
        r_next[i, :] = 2*r_current[i, :] - r_previous[i, :] + F[i, :]*dt**2
        del_r[i, :] = r_next[i, :] - r_previous[i, :]

        if any(r_current[i,:] != r_current_pbc[i,:]):
            # particle i went across the boundary in the previous time step
            r_current[i,:] = r_current_pbc[i,:]
            r_next[i,:] = r_next[i,:]%L

    v_current = del_r/(2*dt)

    return r_current, v_current, r_next, F



"""
Parameters
"""
N = 30      # particle number
L = 10.0    # box length
rc = 2.5    # cutoff-length
T = 1000 # simulation steps
dt = 10**-3 # time step

np.random.seed(42)
"""
Initialization
"""
energy_arr = np.zeros(T)

r_current = np.zeros((N, 3)) # particle positions at time t0
for i in range(N):
    r_current[i] = np.random.rand(3)*L*0.8

v_current = np.random.normal(0, 5, (N, 3)) # particle velocities at time t0
r_next = r_current + v_current*dt # particle positions at time t0+dt

vtk_writer = vtktools.VTK_XML_Serial_Unstructured()
for t in tqdm(range(T)):
    r_current, v_current, r_next, F_ij = stepVerlet(r_current, r_next)
    
    energy_arr[t] = energy(r_next%L, v_current)

    r_current = r_current%L
    r_x = r_current[:, 0]
    r_y = r_current[:, 1]
    r_z = r_current[:, 2]
    F_x = F_ij[:, 0]
    F_y = F_ij[:, 1]
    F_z = F_ij[:, 2]
    vtk_writer.snapshot(os.path.join(data_dir, "MD"+str(t)+".vtu"), r_x, r_y, r_z, \
        x_force=F_x, y_force=F_y, z_force=F_z)

vtk_writer.writePVD(os.path.join(data_dir, "MD.pvd"))

print(F_ij)

"""
Plotting the system energy
"""
plt.figure()
plt.plot(energy_arr)
plt.ylim(0, 1.1*np.max(energy_arr))
plt.xlabel('Timesteps')
plt.ylabel('Energy')
plt.savefig('energy.png')
plt.show()
