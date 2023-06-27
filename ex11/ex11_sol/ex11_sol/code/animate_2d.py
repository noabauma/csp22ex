import argparse, math, random
import matplotlib.pyplot as plt
from src.edmd_2d import edmd_2d


def direct_disks_box(N: int, sigma: float):
    """Generate a valid random initial condition for hard-disks
    by direct sampling.
    
    Keyword arguments:
    N -- number of disks
    sigma -- disk radius

    Returns:
    pos (list of list of floats) -- position vectors of all disks
    """

    condition = False
    while condition == False:
        pos = [[random.uniform(sigma, 1.0 - sigma),
              random.uniform(sigma, 1.0 - sigma)]]
        for k in range(1, N):
            a = [random.uniform(sigma, 1.0 - sigma),
                 random.uniform(sigma, 1.0 - sigma)]
            min_dist = min(
                math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) for b in pos)
            if min_dist < 2.0 * sigma:
                condition = False
                break
            else:
                pos.append(a)
                condition = True
    return pos


def snapshot(t: float, pos: list, vel: list, colors: list, arrow_scale: float=.05):
    """Animate the 2-dimensional simulation."""

    plt.cla()
    plt.axis([0, 1, 0, 1])
    plt.setp(plt.gca(), xticks=[0, 1], yticks=[0, 1])
    for (x, y), (dx, dy), c in zip(pos, vel, colors):
        dx *= arrow_scale
        dy *= arrow_scale
        circle = plt.Circle((x, y), radius=sigma, fc=c)
        plt.gca().add_patch(circle)
        plt.arrow(x, y, dx, dy, fc="k", ec="k",
                  head_width=arrow_scale/2, head_length=arrow_scale/2)
    plt.text(.5, 1.03, 't = %.2f' % t, ha='center')
    plt.pause(0.00001)

    if t == 0:
        plt.savefig("2-d_snapshot.pdf")

plt.subplots_adjust(left=0.10, right=0.90, top=0.90, bottom=0.10)
plt.gcf().set_size_inches(6, 6)
img = 0

"""Run the simulation for specific initial conditions of four hard-disks."""

parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, default=40)
parser.add_argument('-eps', type=float, default=1.)
parser.add_argument('-L', type=float, default=400.)
parser.add_argument('-r', type=float, default=0.03)

args = parser.parse_args()
N = args.N
sigma = args.r
restitution_coefficient = args.eps

# Set the initial conditions for all particles (and colors for visualisation)
pos = direct_disks_box(N, sigma)
vel = [[random.gauss(0, 1) for _ in range(2)] for _ in range(N)]
colors = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])
          for i in range(N)]

edmd = edmd_2d(pos, vel, N, sigma=sigma)

t = 0.0
dt = 0.01     # dt=0 corresponds to event-to-event animation
n_steps = 600

next_event, next_event_arg = edmd.compute_next_event()
snapshot(t, pos, vel, colors)
for step in range(n_steps):
    plt.clf()
    if dt:
        next_t = t + dt
    else:
        next_t = t + next_event
    while t + next_event <= next_t:
        t += next_event
        for k, i in edmd.singles: 
            pos[k][i] += vel[k][i] * next_event
        edmd.compute_new_velocities(next_event_arg, eps_n=restitution_coefficient)
        next_event, next_event_arg = edmd.compute_next_event()
    remain_t = next_t - t
    for k, i in edmd.singles: 
        pos[k][i] += vel[k][i] * remain_t
    t += remain_t
    next_event -= remain_t
    snapshot(t, pos, vel, colors)
