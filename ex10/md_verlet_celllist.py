import os

from requests import delete
import vtktools
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from tqdm import tqdm

np.random.seed(42)

data_dir = "simu"
if not os.path.exists(data_dir):
    os.makedirs(data_dir)

"""
Global Parameters
"""
N = 64      # particle number
L = 10.0    # box length
rc = 2.5    # cutoff-length
rc2 = rc*rc # cutoff squared
n_timesteps = 1000   # simulation steps
dt = 10**-3 # time step
m = 1.0     # particle mass [kg]
dt_2  = dt/2
dt2_2 = dt*dt/2

k_B = 1.0 # Boltzmann constant [m^2 kg s^-2 K^-1]
T   = 100.0 # system temperature [K]

#LJ parameters
epsilon = 1.0
sigma = 1.0
sigma2 = sigma*sigma
V_c = 4*((1/rc)**12-(1/rc)**6) #cutoff potential

#Nosé-Hoover Thermostat parameters
Q = 1000
dt_2Q = dt_2/Q
N3_1_kBT = (3*N+1)*k_B*T
kB_3N_1_inv = 1.0/(k_B * 3*(N-1))


#celllist properties
l = int(L/rc)
l2 = l*l
n_cells = l2*l
cellsize = L/l

def update_celllist(r_current):

    celllist = [ [] for _ in range(n_cells) ]   #empty celllist

    for idx in range(N):
        r = r_current[idx,:]
        i_x = int(r[0]/cellsize)
        i_y = int(r[1]/cellsize)
        i_z = int(r[2]/cellsize)

        celllist[i_x*l2 + i_y*l + i_z].append(idx)

    return celllist

def get_cell_coord(cell_idx):
    cell_idx_yz = cell_idx % l2

    x = int(cell_idx / l2)
    y = int(cell_idx_yz / l)
    z = cell_idx_yz % l

    return x, y, z

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

def LJ_force(r_vec, r2):
    s_r_2 = (sigma2)/r2

    s_r_4 = s_r_2*s_r_2

    s_r_6 = s_r_4*s_r_2

    s_r_8 = s_r_4*s_r_4

    s_r_12 = s_r_8*s_r_4

    return 24.0*epsilon*(2.0*s_r_12 - s_r_6)*r_vec/r2


def LJ_force_potential(r_vec, r2):
    s_r_2 = (sigma2)/r2

    s_r_4 = s_r_2*s_r_2

    s_r_6 = s_r_4*s_r_2

    s_r_8 = s_r_4*s_r_4

    s_r_12 = s_r_8*s_r_4

    potential = 4.0*epsilon*(s_r_12 - s_r_6)

    force = 24.0*epsilon*(2.0*s_r_12 - s_r_6)*r_vec/r2

    return force, potential


def initial_force(r_current):
    celllist = update_celllist(r_current)

    F = np.zeros((N,3))
    for cell_idx in range(n_cells):
        main_cell = celllist[cell_idx]   #main cell
        
        cell_idx_x, cell_idx_y, cell_idx_z  = get_cell_coord(cell_idx)

        for idx_i in range(len(main_cell)):
            i = main_cell[idx_i]
            for idx_x in [(cell_idx_x-1+l)%l, cell_idx_x, (cell_idx_x+1)%l]:
                for idx_y in [(cell_idx_y-1+l)%l, cell_idx_y, (cell_idx_y+1)%l]:
                    for idx_z in [(cell_idx_z-1+l)%l, cell_idx_z, (cell_idx_z+1)%l]:
                        other_cell = celllist[idx_x*l2 + idx_y*l + idx_z]
                        for idx_j in range(len(other_cell)):
                            j = other_cell[idx_j]

                            if i < j:
                                r_vec = r_rel_pbc(r_current[i,:], r_current[j,:])
                                r2 = r_vec.dot(r_vec)
                                if r2 < rc2:
                                    force = LJ_force(r_vec, r2)

                                    F[i,:] += force
                                    F[j,:] -= force

    return F


def step(xi_current, v_current, r_current, F):
    """
    Args:
        r_previous: Particle positions at time t-dt
        r_current: Particle positions at time t

    Returns:
        Updated positions as well as velocities and forces according to the
        Verlet scheme
    """

    E_kin = m*np.sum(np.square(v_current))
            

    #velocity verlet with Nosé-Hoover thermostat

    xi_v = np.multiply(xi_current, v_current.T).T
    
    #full step
    r_next = r_current + v_current*dt + (F/m + xi_v)*dt2_2

    r_next %= L

    celllist = update_celllist(r_next)

    #half step
    v_next  = v_current  + (F/m - xi_v)*dt_2
    xi_next = xi_current + (E_kin - N3_1_kBT)*dt_2Q

    E_kin = m*np.sum(np.square(v_next))
    E_pot = 0.0

    F = np.zeros((N,3))
    for cell_idx in range(n_cells):
        main_cell = celllist[cell_idx]   #main cell
        
        cell_idx_x, cell_idx_y, cell_idx_z  = get_cell_coord(cell_idx)

        for idx_i in range(len(main_cell)):
            i = main_cell[idx_i]
            for idx_x in [(cell_idx_x-1+l)%l, cell_idx_x, (cell_idx_x+1)%l]:
                for idx_y in [(cell_idx_y-1+l)%l, cell_idx_y, (cell_idx_y+1)%l]:
                    for idx_z in [(cell_idx_z-1+l)%l, cell_idx_z, (cell_idx_z+1)%l]:
                        other_cell = celllist[idx_x*l2 + idx_y*l + idx_z]
                        for idx_j in range(len(other_cell)):
                            j = other_cell[idx_j]

                            if i < j:
                                r_vec = r_rel_pbc(r_next[i,:], r_next[j,:])
                                r2 = r_vec.dot(r_vec)
                                if r2 < rc2:
                                    force, potential = LJ_force_potential(r_vec, r2)

                                    E_pot += potential - V_c
                                    F[i,:] += force
                                    F[j,:] -= force

    xi_next += (E_kin - N3_1_kBT)*dt_2Q

    #v_next = (v_next + dt_2*F/m)/(1 + dt_2*xi_next)
    v_next = np.multiply(1.0/(1.0 + dt_2*xi_next),(v_next + dt_2*F/m).T).T


    return xi_next, v_next, r_next, F, 0.5*E_kin + E_pot, kB_3N_1_inv*E_kin





def main():
    """
    Initialization
    """
    energy_arr = np.zeros(n_timesteps)
    T_inst_arr = np.zeros(n_timesteps)

    xi_current = np.zeros(N)

    r_current = L*np.random.rand(N,3)


    # https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
    var = np.sqrt(k_B*T/m)*0.8
    v_current = np.random.normal(scale=var, size=(N,3))


    F_ij = initial_force(r_current)

    # Run the time evolution for `n_timesteps` steps:
    vtk_writer = vtktools.VTK_XML_Serial_Unstructured()
    for t in tqdm(range(n_timesteps)):

        xi_current, v_current, r_current, F_ij, energy_arr[t], T_inst_arr[t] = step(xi_current, v_current, r_current, F_ij)

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
    fig, ax = plt.subplots(2)

    ax[0].plot(np.linspace(0,n_timesteps,2), T*np.ones(2), 'r--')
    ax[0].plot(T_inst_arr)
    ax[0].set_ylim(0, 1.1*np.max(T_inst_arr))
    ax[0].set_ylabel('Temperature')

    ax[1].plot(energy_arr)
    ax[1].set_ylim(0, 1.1*np.max(energy_arr))
    ax[1].set_xlabel('Timesteps')
    ax[1].set_ylabel('Energy')

    ax[0].set_title('T = ' + str(T) + ', Q = ' + str(Q))
    plt.savefig('temperature&energy.png')
    plt.show()


if __name__ == "__main__":
    main()