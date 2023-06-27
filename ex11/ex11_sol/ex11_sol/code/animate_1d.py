import random, argparse
import numpy as np
import matplotlib.pyplot as plt
from src.edmd_1d import *


def ed_animation(n_steps: int, R: float, N: int, L: float, 
        restitution_coefficient: float, dt: float=0.001, verbose: bool=True):
    #c_data = np.random.rand(N)
    c_data = np.linspace(0,1,N)

    x, v = generate_particles(N, R, L, d=1)  # randomly initialise particles
    t = 0  # initial time (time since last collision)

    t_c, i_c, wall = collision_finder_wall(x, v, N, L, R=R, verbose=verbose)
    for step in range(n_steps):
        if dt:
            next_t = t + dt
        else:
            next_t = t + t_c

        while t + t_c <= next_t:
            t += t_c
            x = position_update(x, v, t_c, N)  # update the positions
            v = collision_velocity_update(
                v, i_c, N, wall=wall, eps_n=restitution_coefficient)

            t_c, i_c, wall = collision_finder_wall(x, v, N, L, R=R, verbose=verbose)

        remain_t = next_t - t
        x = position_update(x, v, remain_t, N)  # update the positions

        t += remain_t
        t_c -= remain_t

        if step % 20 == 0:
            snapshot_beads(x, N, sigma, c_data, L=L)


"""Run simulation for specific initial conditions of four hard-disks."""


parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, default=15)
parser.add_argument('-eps', type=float, default=1.)
parser.add_argument('-L', type=float, default=400.)
parser.add_argument('-r', type=float, default=1.)

args = parser.parse_args()

sigma = args.r
L = args.L
N = args.N
restitution_coefficient = args.eps

dt = 0.05
n_steps = int(10e6)

fig = plt.figure(figsize=(L, 1.5))
a = ed_animation(n_steps, sigma, N, L,
                 restitution_coefficient, dt=dt)
