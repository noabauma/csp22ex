import math

def min_arg(l: list):
    """Determines the value and the location of the
        minimal element in a list.

    Keyword argument:
    l -- list to be searched

    Returns:
    min(l) (float) -- miniemal element
    argmin(l) (integer) -- location of the minimal element
    """

    return min(zip(l, range(len(l))))


class edmd_2d():
    def __init__(self, pos: list, vel: list, N: int, sigma: float=0.15, L: float=1.0):
        self.sigma = sigma # disk radius
        self.L = L

        self.pos = pos
        self.vel = vel

        self.N = N # number of disks
        # List tuples of disk ids and directions:
        self.singles = [(k, i) for k in range(N) for i in range(2)]
        # List tuples of ids representing pairs of disks:
        self.pairs = [(k, l) for k in range(N) for l in range(N) if k < l]


    def wall_time(self, pos_i: float, vel_i: float):
        """Computes the collision time with a boundary wall
        in the i'th direction.

        Keyword arguments:
        pos_i -- i'th position component for the particle
        vel_i -- i'th velocity component for the particle

        Returns:
        del_t (float) -- collsion time with a wall in i'th direction
        """

        if vel_i > 0.0:
            del_t = (self.L - self.sigma - pos_i) / vel_i
        elif vel_i < 0.0:
            del_t = (pos_i - self.sigma) / abs(vel_i)
        else:
            del_t = float('inf')

        return del_t


    def pair_time(self, pos_a: list, vel_a: list, pos_b: list, vel_b: list):
        """Computes the collision time of a pair of disks a and b.

        Keyword arguments:
        pos_a(b) (list of floats) -- position vector for particle a (b)
        vel_a(b) (list of floats) -- velocity vector for particle a (b)

        Returns:
        del_t (float) -- the smallest time until collsion
        """

        del_x = [pos_b[0] - pos_a[0], pos_b[1] - pos_a[1]]
        del_x_sq = del_x[0] ** 2 + del_x[1] ** 2

        del_v = [vel_b[0] - vel_a[0], vel_b[1] - vel_a[1]]
        del_v_sq = del_v[0] ** 2 + del_v[1] ** 2

        scal = del_v[0] * del_x[0] + del_v[1] * del_x[1]
        Upsilon = scal ** 2 - del_v_sq * (del_x_sq - 4.0 * self.sigma ** 2)
        if Upsilon > 0.0 and scal < 0.0:
            del_t = - (scal + math.sqrt(Upsilon)) / del_v_sq
        else:
            del_t = float('inf')

        return del_t


    def compute_next_event(self):
        """Determines the earliest collision event.

        Returns:
        collision_time (float) -- collision time for the earliest event
                                    and the index for the next event
        collision_event (integer) -- index specifying the collision event:
                                        either the disk-wall collision or 
                                        the specific disk-disk collision.

        Example:
        if collision_event < len(singles): 
            # the earliest collision is with a wall
        else:
            # the earliest collision is between pairs of disks
        """

        wall_times = [self.wall_time(self.pos[k][i], self.vel[k][i]) 
                        for k, i in self.singles]
        pair_times = [self.pair_time(self.pos[k], self.vel[k], self.pos[l], self.vel[l]) 
                        for k, l in self.pairs]

        collision_time, collision_event = min_arg(wall_times + pair_times)
        return collision_time, collision_event


    def compute_new_velocities(self, next_event_arg: tuple, eps_n: float=1):
        """Updates the postition and velocity vectors after a collision event

        Keyword arguments:
        next_event_arg (tuple of a float and an integer) -- (collision time, event)
        eps_n -- restitution coefficient
        """

        if next_event_arg < len(self.singles): # if the earliest collision is with a wall
            collision_disk, direction = self.singles[next_event_arg]
            self.vel[collision_disk][direction] *= -1.0
        else: # if the earliest collision is between pairs of disks
            a, b = self.pairs[next_event_arg - len(self.singles)]
            del_x = [self.pos[b][0] - self.pos[a][0], self.pos[b][1] - self.pos[a][1]]
            abs_x = math.sqrt(del_x[0] ** 2 + del_x[1] ** 2)
            e_perp = [c / abs_x for c in del_x]
            del_v = [self.vel[b][0] - self.vel[a][0], self.vel[b][1] - self.vel[a][1]]
            scal = del_v[0] * e_perp[0] + del_v[1] * e_perp[1]

            for k in range(2):
                self.vel[a][k] += e_perp[k] * scal
                self.vel[b][k] -= e_perp[k] * scal

                # dissipate energy due to resitution
                self.vel[a][k] *= eps_n
                self.vel[b][k] *= eps_n
