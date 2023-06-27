import os
import vtktools
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from tqdm import tqdm


data_dir = "simu"
if not os.path.exists(data_dir):
    os.makedirs(data_dir)

"""
Global Parameters
"""
N = 30      # particle number
L = 10.0    # box length
rc = 2.5    # cutoff-length
rc2 = rc*rc # cutoff squared
T = 1000   # simulation steps
dt = 10**-3 # time step
m = 1.0     # particle mass [kg]
dt2_m = dt*dt/m # dt squared divided by mass
dt_m = dt/m
dt2_inv = 1.0/(2.0*dt)

np.random.seed(42)

#LJ parameters
epsilon = 1.0
sigma = 1.0

V_c = 4*((1/rc)**12-(1/rc)**6) #cutoff potential


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
        if abs(r_k) > L*0.5:
            r_vec[k] = -np.sign(r_k)*(L - abs(r_k))
    
    return r_vec
    

def force_potential(r_vec, r2):
    s_r_2 = (sigma*sigma)/r2

    s_r_4 = s_r_2*s_r_2

    s_r_6 = s_r_4*s_r_2

    s_r_8 = s_r_4*s_r_4

    s_r_12 = s_r_8*s_r_4

    potential = 4.0*epsilon*(s_r_12 - s_r_6)

    force = 24.0*epsilon*(2.0*s_r_12 - s_r_6)*r_vec/r2

    return force, potential





def stepVerlet(r_previous, v_current, r_current):
    """
    Args:
        r_previous: Particle positions at time t-dt
        r_current: Particle positions at time t

    Returns:
        Updated positions as well as velocities and forces according to the
        Verlet scheme
    """

    E_kin = 0.5*m*np.sum(np.square(v_current))
    E_pot = 0.0
    
    F = np.zeros((N, 3))
    # TODO: compute the total force (=acceleration) acting on each particle
    for i in range(N):
        for j in range(i+1,N):
            r_vec = r_rel_pbc(r_current[i,:], r_current[j,:])
            r2 = r_vec.dot(r_vec)
            if r2 < rc2:
                force, potential = force_potential(r_vec, r2)

                E_pot += potential - V_c

                F[i,:] += force
                F[j,:] -= force
    

    #verlet step
    """
    r_next = 2.0*r_current - r_previous + dt2_m*F

    v_current = (r_next - r_previous)*dt2_inv

    for i in range(N):
        if any((r_current[i,:] > L) & (r_next[i,:] > L)):
        #if any(r_current[i,:] > L):
                r_current[i,:] %= L
                r_next[i,:] %= L
    """

    #velocity verlet
    v_current += dt_m*F
    r_next = r_current + v_current*dt

    r_next %= L

    return r_current, v_current, r_next, F, E_kin + E_pot





def main():
    """
    Initialization
    """
    energy_arr = np.zeros(T)

    r_current = L*np.random.rand(N,3)*0.8

    #this is normally how you initiallize the velocity of the particles (it is temperature dependent!):
    # https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
    k_b  = 1.0 # Boltzmann constant [m^2 kg s^-2 K^-1]
    temp = 25.0 # system temperature [K]
    var = np.sqrt(k_b*temp/m)

    v_current = np.random.normal(scale=var, size=(N,3))

    r_next = r_current + v_current*dt # particle positions at time t0+dt

    # Run the time evolution for `T` steps:
    vtk_writer = vtktools.VTK_XML_Serial_Unstructured()
    for t in tqdm(range(T)):
        r_current, v_current, r_next, F_ij, energy_arr[t] = stepVerlet(r_current, v_current, r_next)

        #r_current = r_current%L
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
    #print(energy_arr)
    plt.figure()
    plt.plot(energy_arr)
    plt.ylim(0, 1.1*np.max(energy_arr))
    #plt.ylim(0,1100)
    plt.xlabel('Timesteps')
    plt.ylabel('Energy')
    plt.savefig('energy2.png')
    plt.show()


if __name__ == "__main__":
    main()