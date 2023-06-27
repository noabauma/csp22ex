## Event-driven molecular dynamics
This repositiory contains the code for the event-driven molecular dynamics simulations of isotropic particles with hard-core interactions.

* `src/edmd1d.py` contains the source code for 1-dimensional rods. It can be used by running `animate1d.py` from the terminal with arguments 

  * `-N` for the number of beads
  * `-L` for the size of the box
  * `-r` for the width of the rods
  * `-eps` for the restitution coefficient (`=1`: no energy dissipation, `=0`: completely inelastic collisions)

  Also, `restitution.ipynb` notebook generates the "phase diagram" by plotting the total effective restitution coefficient of the rod system as a function of the restitution coefficient of collisions.

* `src/edmd2d.py` is the source code for 2-dimensional hard-disks. Similarly it can be used by running `animate2d.py` with the same arguments, but now `-r` stands for the radius of the disks. 

Note that although the animations are generated for closed boundary conditions (with perfectly inelastic collisions with walls), the source codes also include the functions for determining collision events and performing velocity updates with periodic boundary conditions.