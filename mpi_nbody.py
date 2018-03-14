# ====================================================================
# Author 				: swc21
# Date 					: 2018-03-14 09:42:27
# Project 				: ClusterFiles
# File Name 			: mpi_nbody
# Last Modified by 		: swc21
# Last Modified time 	: 2018-03-14 10:28:53
# ====================================================================
#
#<Begin: Imports>
from mpi4py import MPI
import numpy as np
import sys
#<End: Imports>
#<Begin: MPI Init>
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
TIMESTEPS = 50
if len(sys.argv) > 1:
    try:
        new_timesteps = int(sys.argv[1])
        assert isinstance(new_timesteps, int)
        assert new_timesteps > 1
        TIMESTEPS = new_timesteps
    except Exception:
        pass
#<End: MPI Init>
#<Begin: Particle Init>


def random_particle(rank_seed=rank):
    x_y_z_m = np.random.rand(4)
    x_y_z_m[:3] *= 100.0
    vx_vy_vz = (np.random.rand(3)-0.5)*2.0
    ax_ay_az = np.zeros(3)
    return x_y_z_m, vx_vy_vz, ax_ay_az


my_xyzm, my_v, my_a = random_particle()
#<End: Particle Init>
#<Begin: Physics>


# =my_xyzm):
def change_in_a_from_xyzm(other_xyzm, xyzm, G=1e1, tolerance=1e1):
    # TODO
    if np.all(other_xyzm == xyzm):
        return np.zeros(3)
    displacement = other_xyzm[:3]-xyzm[:3]
    distance = np.sqrt(np.square(displacement).sum())
    if distance < tolerance:
        return np.zeros(3)
    return G*other_xyzm[3]*displacement/distance**3


def update_state(many_other_xyzm, xyzm=my_xyzm, v=my_v, a=my_a, DT=1.0):
    a = np.zeros(3)
    # TODO May Want To Do This Destructively
    for other_xyzm in many_other_xyzm:
        a += change_in_a_from_xyzm(other_xyzm, xyzm)
    v += a*DT
    xyzm[:3] += v*DT
    return xyzm, v, a  # TODO May Not Be Needed


#<End: Physics>
#<Begin: Timestep>
if not rank:
    STATES_TO_SAVE = 300
    assert 1000 >= STATES_TO_SAVE >= 2
    STATES = []
for TIMESTEP in range(TIMESTEPS):
    # TODO Calculate CM for my internal points and do internal nbody
    if not rank:  # Test Function
        if not TIMESTEP:
            print 'INITIAL:'
        print 'my_xyzm', my_xyzm
    all_xyzm = comm.allgather(my_xyzm)  # Warning: This Contains Self
    if (not rank) and (not TIMESTEP % (TIMESTEPS//STATES_TO_SAVE)):  # Test Function
        # print 'all_xyzm', all_xyzm
        # while len(STATES) >= STATES_TO_SAVE:
        # STATES.pop()
        STATES.append(all_xyzm)
    my_xyzm, my_v, my_a = update_state(all_xyzm)
#<End: Timestep>
if not rank:  # Test Function
    print 'FINAL:'
    print 'my_xyzm', my_xyzm
    # print 'all_xyzm', all_xyzm
    import os
    np.save(os.path.join('SHARED', 'mpi_nbody_out.npy'), np.array(STATES))
