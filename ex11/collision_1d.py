import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from tqdm import tqdm
import sys

np.random.seed(38)

#global parameters
N = 40 # number of particles
T = 1000 # number of timesteps
L = 100 # system size [0, L]

m = 1.0 # mass of the particle

save_every = 40 # only save every so frame

n_frames = 1000

dt = 0.01 #timestep

R = 1    # radius of the particles
R2 = (2*R)**2 #
R_wall = R**2

e = 0.9 # restitution coefficient

#this implementation only works for 1d case!
def global_collision_time(r, v):

    # particles on the left boundary
    r_ij = -r[0]
    v_ij = -v[0]

    a = v_ij*v_ij
    b = 2*r_ij*v_ij
    c = r_ij*r_ij - R_wall

    temp1 = (-b + np.sqrt(b*b - 4*a*c))/(2*a)
    temp2 = (-b - np.sqrt(b*b - 4*a*c))/(2*a)

    t_ = 0.0
    if (temp1 > 0) and (temp2 > 0):
        t_ = temp1 if temp1 < temp2 else temp2
    elif (temp1 > 0):
        t_ = temp1
    elif (temp2 > 0):
        t_ = temp2
    else:
        t_ = np.inf

    t_C = [t_, -1]   # collision time, particle i

    for i in range(N-1):        

        #right neighbor
        r_ij = r[i] - r[i+1]
        v_ij = v[i] - v[i+1]

        a = v_ij*v_ij
        b = 2*r_ij*v_ij
        c = r_ij*r_ij - R2

        temp1 = (-b + np.sqrt(b*b - 4*a*c))/(2*a)
        temp2 = (-b - np.sqrt(b*b - 4*a*c))/(2*a)

        t_ = 0.0
        if (temp1 > 0) and (temp2 > 0):
            t_ = temp1 if temp1 < temp2 else temp2
        elif (temp1 > 0):
            t_ = temp1
        elif (temp2 > 0):
            t_ = temp2
        else:
            t_ = np.inf

        if t_ < t_C[0]:
            t_C = [t_, i]

    # particles on the right boundary
    r_ij = r[N-1] - L
    v_ij = v[N-1]

    a = v_ij*v_ij
    b = 2*r_ij*v_ij
    c = r_ij*r_ij - R_wall

    temp1 = (-b + np.sqrt(b*b - 4*a*c))/(2*a)
    temp2 = (-b - np.sqrt(b*b - 4*a*c))/(2*a)

    t_ = 0.0
    if (temp1 > 0) and (temp2 > 0):
        t_ = temp1 if temp1 < temp2 else temp2
    elif (temp1 > 0):
        t_ = temp1
    elif (temp2 > 0):
        t_ = temp2
    else:
        t_ = np.inf

    if t_ < t_C[0]:
        t_C = [t_, N-1]

    return t_C



def main():
    r = np.linspace(3, 97, N)     # initial coordinates
    v = np.random.rand(N)*3 - 1 # initial velocity

    r_tot = np.empty((n_frames,N))
    v_tot = np.empty((n_frames,N))

    r_tot[0,:] = r
    v_tot[0,:] = v

    E_tot = []
    E_tot.append(np.sum(v*v))

    frame = 1
    while frame < n_frames:
        t_C = global_collision_time(r, v)

        # take so many steps until next collision accures
        for step in range(int(t_C[0]/dt)):
            r += v*dt

            if (step % save_every) == 0:
                r_tot[frame,:] = r
                v_tot[frame,:] = v

                frame += 1

                if (frame+1) >= n_frames:
                    frame += 1
                    break

        if (frame+1) >= n_frames:
            break

        # update velocities of the meant collision
        
        i = t_C[1]
        if i == -1:  # check if it is a collision with left wall
            tmp = v[i+1]
            v[i+1] *= -e
            E_tot.append(E_tot[-1] - tmp**2 + v[i+1]**2)
        elif i == N-1: # check if it is a collision with right wall
            tmp = v[i]
            v[i] *= -e
            E_tot.append(E_tot[-1] - tmp**2 + v[i]**2)
        else:             # not a wall
            tmp_i = v[i]
            tmp_j = v[i+1]

            r_ij = r[i] - r[i+1]
            n = r_ij/abs(r_ij)

            v_ij = (v[i] - v[i+1])*n

            v[i]   -= v_ij*n
            v[i+1] += v_ij*n

            v[i]   *= e
            v[i+1] *= e

            E_tot.append(E_tot[-1] - tmp_i**2 + v[i]**2 - tmp_j**2 + v[i+1]**2)
    

    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes(xlim=(0, L), ylim=(-1,1))
    line, = ax.plot([], [], 'o')

    # initialization function: plot the background of each frame
    def init():
        line.set_data([], [])
        return line,

    # animation function.  This is called sequentially
    def animate(i):
        x = r_tot[i,:]
        y = np.zeros(N)
        line.set_data(x, y)
        return line,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=n_frames, interval=5, blit=True)
    
    #anim.save('animation_1d.mp4', fps=120)

    plt.show()


    E_tot = np.array(E_tot)
    E_tot *= 0.5*m

    print('effective restitution coefficient e_eff = ', np.sqrt(E_tot[-1]/E_tot[0]))

    len_E = E_tot.shape[0]
    plt.plot(np.arange(len_E), np.round(E_tot, 2))
    plt.xlabel('timestep')
    plt.ylabel('Energy')
    plt.title('total Energy for e = ' + str(e))
    plt.show()




if __name__ == "__main__":
    main()