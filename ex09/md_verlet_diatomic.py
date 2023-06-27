import os

from requests import delete
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
n = 6       # 
N_2 = n**3  # half N 
N = N_2*2   # particle number
L = 10.0    # box length
rc = 2.5    # cutoff-length
rc2 = rc*rc # cutoff squared
T = 1000    # simulation steps
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

#lagrange
d = epsilon/10.0    #the forced distance constrain
d2 = d*d
dt2_2_inv = dt2_inv/dt


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
    

def LJ_force_potential(r_vec, r2):
    s_r_2 = (sigma*sigma)/r2

    s_r_4 = s_r_2*s_r_2

    s_r_6 = s_r_4*s_r_2

    s_r_8 = s_r_4*s_r_4

    s_r_12 = s_r_8*s_r_4

    potential = 4.0*epsilon*(s_r_12 - s_r_6)

    force = 24.0*epsilon*(2.0*s_r_12 - s_r_6)*r_vec/r2

    return force, potential


def stepVerlet(r_previous, v_current, r_current, celllist):
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

                            #if i < j+1:
                            if i < j-1:
                                r_vec = r_rel_pbc(r_current[i,:], r_current[j,:])
                                r2 = r_vec.dot(r_vec)
                                if r2 < rc2:
                                    force, potential = LJ_force_potential(r_vec, r2)

                                    E_pot += potential - V_c
                                    F[i,:] += force
                                    F[j,:] -= force
            
    
    
    #verlet step
    r_next = 2.0*r_current - r_previous + dt2_m*F

    lambda_ = np.empty(N)
    for i in range(N_2):
        r_12 = r_current[2*i,:] - r_current[2*i+1,:]
        r_12_2 = r_12.dot(r_12)

        r_12_tilde = r_next[2*i,:] - r_next[2*i+1,:]
        r_12_2_tilde = r_12_tilde.dot(r_12_tilde)

        dt_temp = dt2_2_inv/r_12_2

        r_r_tilde = r_12.dot(r_12_tilde)
        #temp = np.sqrt(r_r_tilde**2 - r_12_2*(r_12_2_tilde - d2))
        temp = np.sqrt(r_r_tilde**2 - d2*(r_12_2_tilde - d2))      #should be the same!

        lambda_plus  = dt_temp*(-r_r_tilde + temp)
        lambda_minus = dt_temp*(-r_r_tilde - temp)

        lambda_magn = lambda_plus if abs(lambda_plus) < abs(lambda_minus) else lambda_minus
        
        lambda_[2*i]   = lambda_magn
        lambda_[2*i+1] = lambda_magn

    g_ = r_next[::2, :] - r_next[1::2, :]
    g  = np.repeat(g_, 2, axis=0)
    g[1::2, :] *= -1

    #print(r_next.shape, lambda_.shape, dt2_m, g.shape)

    for i in range(N):  #TODO: this can be transformed into numpy operators
        r_next[i,:] += lambda_[i]*dt2_m*g[i,:]

    #r_next += lambda_*dt2_m*g

    v_current = (r_next - r_previous)*dt2_inv

    for i in range(N):
        if any((r_current[i,:] > L) | (r_next[i,:] > L)):
        #if any(r_current[i,:] > L):
                r_current[i,:] %= L
                r_next[i,:] %= L

    return r_current, v_current, r_next, F, E_kin + E_pot





def main():
    """
    Initialization
    """
    energy_arr = np.zeros(T)

    #shifting the paired particle by vectors of length d
    d_vec = np.random.rand(N_2,3)
    d_norm = np.linalg.norm(d_vec, axis=1)
    for i in range(N_2): #TODO: this can be transformed into numpy operators
        d_vec[i,:] /= d_norm[i]
    d_vec *= 1.5*d


    r_current = np.empty((N,3))
    displacement = L/float(n+1)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                idx = (i*n**2 + j*n + k)
                idx2 = 2*idx
                r_current[idx2,:] = np.array([(i+1)*displacement, (j+1)*displacement, (k+1)*displacement])
                r_current[idx2+1,:] = r_current[idx2,:] + d_vec[idx,:]

    r_current %= L

    #r_current = L*np.random.rand(N,3)


    #this is normally how you initiallize the velocity of the particles (it is temperature dependent!):
    # https://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution
    k_b  = 1.0 # Boltzmann constant [m^2 kg s^-2 K^-1]
    temp = 2.0 # system temperature [K]
    var = np.sqrt(k_b*temp/m)

    v_current = np.random.normal(scale=var, size=(N,3))

    r_next = r_current + v_current*dt # particle positions at time t0+dt

    # Run the time evolution for `T` steps:
    vtk_writer = vtktools.VTK_XML_Serial_Unstructured()
    for t in tqdm(range(T)):

        celllist = update_celllist(r_next)

        r_current, v_current, r_next, F_ij, energy_arr[t] = stepVerlet(r_current, v_current, r_next, celllist)

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
    plt.savefig('energy.png')
    plt.show()


if __name__ == "__main__":
    main()