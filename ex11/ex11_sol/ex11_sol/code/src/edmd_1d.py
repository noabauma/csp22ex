import random
from tqdm.notebook import tqdm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
rc('font', size=10)
rc('axes', labelsize=10)


def dict_to_array(d: dict):
    # array storing vectors for all N particles
    d_array = np.array([v for v in d.values()])
    return d_array

def vec_rel(vec, i: int, j: int): return (np.longdouble(vec[i]-vec[j]))


def snapshot(pos: dict, colors, L: float=1.0, sigma: float=0.1):
    # begin of graphical output
    plt.cla()
    bluesquare = plt.Rectangle(
        (sigma, 0), L - 2 * sigma, 0.33 * L, fc='b', alpha=0.2)
    plt.gca().add_patch(bluesquare)
    for pin in pos:
        whiterec = plt.Rectangle((pin - 2 * sigma, 0),
                                 4 * sigma, 0.33 * L, fc='w', ec='w')
        plt.gca().add_patch(whiterec)
    for i, pin in enumerate(pos):
        redrec = plt.Rectangle((pin - sigma, 0), 2 *
                               sigma, 0.33 * L, fc=colors[i])
        plt.gca().add_patch(redrec)
    plt.axis('scaled')
    plt.axis([0, L, 0, 0.2 * L])
    plt.xlabel('$x$', fontsize=14)
    plt.pause(10e-4)


def snapshot_beads(pos: dict, N: int, sigma: float, colors: np.ndarray, L: float=1.0):
    plt.clf()
    plt.scatter(np.reshape(dict_to_array(pos), N),
                np.zeros(N), marker='o', c=cm.jet(colors),
                linewidth=(90*sigma)**2/4/L)

    plt.xlabel('$x$')
    plt.xlim(0, L)
    plt.ylim(-0.02, 0.02)
    plt.yticks([])
    plt.pause(10e-4)

def generate_particles(N: int, R: float, L: float, d):
    coordinates = {}
    velocities = {}  # initialise coordinates and velocities
    coordinates_values = np.array([2*R*i + x for i, x in
        enumerate(sorted(random.sample(np.linspace(R, L-R-(N-1)*2*R, 1000).tolist(), N)))])
    for i in range(N):
        coordinates[i] = np.sort(coordinates_values)[i]
        velocities[i] = np.random.uniform(8, 10, d)
    return coordinates, velocities  # returns all N coordinates and velocities


# determine the time until the next collision and the collision time
def collision_finder_periodic(x: np.ndarray, v: np.ndarray, N: int, 
                                    L: float, R=0, verbose: bool = False):
    # initialise the set of all neighbouring encounter times
    t_ij = np.zeros(N, dtype=np.longdouble)
    for i in range(N):
        j = (i+1) % N  # in 1-d only nearest neighbours can collide
        x_ji = vec_rel(x, j, i) % L
        # relative positions and velocities of neighbouring particles in pbc
        v_ij = vec_rel(v, i, j)
        t_ij[i] = x_ji/v_ij-2*R/abs(v_ij)
        if t_ij[i] < 0:
            t_ij[i] = float('inf')
    t_c = np.min(t_ij)
    i_c = np.argmin(t_ij)
    if verbose: print(t_c, i_c)
    return t_c, i_c


def collision_finder_wall(x: np.ndarray, v: np.ndarray, N: int, 
            L: float, R: float=0, threshold=10e-20, verbose: bool = False):
    # initialise the set of all neighbouring encounter times
    t_ij = np.ones(N, dtype=np.longdouble)*float('inf')
    t_w = np.zeros(N, dtype=np.longdouble)
    wall = False
    for i in range(N):
        if i != N-1:
            j = i+1  # in 1-d only nearest neighbours can collide
            x_ji = vec_rel(x, j, i)
            # relative positions and velocities of neighbouring particles
            v_ij = vec_rel(v, i, j)
            t_ij[i] = x_ji/v_ij-2*R/abs(v_ij)
            if t_ij[i] < 0:
                t_ij[i] = float('inf')

        if v[i] > threshold:
            t_w[i] = -x[i]/v[i]+(L-R)/v[i]
        else:
            t_w[i] = x[i]/abs(v[i])-R/abs(v[i])
    t_c = min(np.min(t_ij), np.min(t_w))
    if np.min(t_w) == t_c:
        wall = True
        i_c = np.argmin(t_w)
    else:
        i_c = np.argmin(t_ij)
    if verbose: print(t_c, i_c)
    return t_c, i_c, wall


def collision_velocity_update(v: np.ndarray, i_c: int, N: int, 
                                wall: bool = False, M=1, eps_n: float = 1):
    i = i_c
    j = (i_c+1) % N
    v_i_old = v[i].copy()
    v_j_old = v[j].copy()
    #update velocities
    if wall == False:
        if np.size(M) == 1:
            v[i] = eps_n*v_j_old
            v[j] = eps_n*v_i_old
        else:
            v[i] = (M[i]-M[j])/(M[i]+M[j])*v_i_old + \
                2*M[j]/(M[i]+M[j])*v_j_old*eps_n
            v[j] = 2*M[i]/(M[i]+M[j])*v_i_old + (M[j]-M[i]) / \
                (M[i]+M[j])*v_j_old*eps_n
    else:
        if np.size(M) == 1:
            v[i] *= -1  # perfect reflection from the walls
    return v


def position_update(x: np.ndarray, v: np.ndarray, t: float, N: int):
    for i in range(N):
        x[i] = (x[i]+t*v[i])
    return x


def ed_simulation(x: np.ndarray, v: np.ndarray, n_steps: int, R: float, N: int, 
            L: float, restitution_coefficient: float = 1, 
            t_max: float = float('inf'), verbose: bool = False):
    t = 0  
    t_c, i_c, wall = collision_finder_wall(x, v, N, L, R=R, verbose=verbose)
    for step in tqdm(range(n_steps), disable=not verbose):
        t += t_c
        x = position_update(x, v, t_c, N)  # update the positions
        v = collision_velocity_update(
            v, i_c, N, wall=wall, eps_n=restitution_coefficient)

        t_c, i_c, wall = collision_finder_wall(x, v, N, L, R=R, verbose=verbose)

        if t > t_max:
            break
    return x, v
    