import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from tqdm import tqdm
import sys

np.random.seed(38)

#global parameters
N_sqrt = 5
N = N_sqrt**2 # number of particles
T = 1000 # number of timesteps
L = 50 # system size [0, L]

save_every = 40 # only save every so frame

n_frames = 1000

dt = 0.01 #timestep

R = 1    # radius of the particles
R2 = (2*R)**2 #
R_wall = R**2

#this implementation only works for 1d case!
def global_collision_time(r, v):
    t_C = [np.inf, np.inf, np.inf]
    for i in range(N):
        for j in range(i+1,N):
            #right neighbor
            r_ij = r[i,:] - r[j,:]
            v_ij = v[i,:] - v[j,:]

            a = v_ij.dot(v_ij)
            b = 2*r_ij.dot(v_ij)
            c = r_ij.dot(r_ij) - R2

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
                t_C = [t_, i, j]

        #the 4 walls separate   top:-1, left:-2, bottom:-3, right:-4

        # top wall
        r_ij = r[i,0] - L
        v_ij = v[i,0]

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
            t_C = [t_, i, -1]

        # left wall
        r_ij = r[i,1] - L
        v_ij = v[i,1]

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
            t_C = [t_, i, -2]

        # bottom wall
        r_ij = r[i,0]
        v_ij = v[i,0]

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
            t_C = [t_, i, -3]

        # right wall
        r_ij = r[i,1]
        v_ij = v[i,1]

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
            t_C = [t_, i, -4]
    
    return t_C


def main():
    r = np.empty((N,2))
    for i in range(N_sqrt):
        for j in range(N_sqrt):
            r[i*N_sqrt + j,:] = [i/N_sqrt*L + 2.5, j/N_sqrt*L + 2.5]
    
    v = np.random.normal(size=(N,2))*0.1 # initial velocity

    r_tot = np.empty((n_frames,N,2))
    v_tot = np.empty((n_frames,N,2))

    r_tot[0,:,:] = r
    v_tot[0,:,:] = v

    frame = 1
    while frame < n_frames:
        t_C = global_collision_time(r, v)

        # take so many steps until next collision accures
        for step in range(int(t_C[0]/dt)):
            r += v*dt

            if (step % save_every) == 0:
                r_tot[frame,:,:] = r
                v_tot[frame,:,:] = v

                frame += 1

                if (frame+1) >= n_frames:
                    frame += 1
                    break

        if (frame+1) >= n_frames:
            break

        # update velocities of the meant collision
        
        i = t_C[1]
        j = t_C[2]
        if (j == -1) or (j == -3):  # check top or bottom wall
            v[i,0] *= -1
        elif (j == -2) or (j == -4): # check left or right wall
            v[i,1] *= -1
        else:             # not a wall
            r_ij = r[i,:] - r[j,:]
            n = r_ij/np.linalg.norm(r_ij)

            v_ij = (v[i] - v[j])*n

            v[i,:] -= v_ij*n
            v[j,:] += v_ij*n    
    

    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = plt.axes(xlim=(0, L), ylim=(0, L))
    line, = ax.plot([], [], 'o', markersize=11)

    # initialization function: plot the background of each frame
    def init():
        line.set_data([], [])
        return line,

    # animation function.  This is called sequentially
    def animate(i):
        x = r_tot[i,:,0]
        y = r_tot[i,:,1]
        line.set_data(x, y)
        return line,

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=n_frames, interval=5, blit=True)
    
    #anim.save('animation_2d.mp4', fps=120)

    plt.show()




if __name__ == "__main__":
    main()