# ====================================================================
# Author 				: swc21
# Date 					: 2018-03-14 09:42:27
# Project 				: ClusterFiles
# File Name 			: mpi_nbody_updated_lowell
# Last Modified by 		: swc21
# Last Modified time 	: 2018-03-14 10:28:54
# ====================================================================
#

#<Begin: Imports>
from mpi4py import MPI
import numpy as np
import sys
#<End: Imports>


#<Begin: Command Line Args>
TIMESTEPS = 100
if len(sys.argv) > 1:
    try:
        new_timesteps = int(sys.argv[1])
        assert isinstance(new_timesteps, int)
        assert 100000 > new_timesteps >= 2, "Safety Measure"
        TIMESTEPS = new_timesteps
    except Exception:
        TIMESTEPS = 100

STATES_TO_SAVE = 10
if len(sys.argv) > 2:
    try:
        new_states_to_save = int(sys.argv[2])
        assert isinstance(new_states_to_save, int)
        assert min(100, TIMESTEPS) >= new_states_to_save >= 2, "Safety Measure"
        STATES_TO_SAVE = new_states_to_save
    except Exception:
        STATES_TO_SAVE = 10

VELOCITY_VERLET = False
if len(sys.argv) > 3:
    try:
        new_velocity_verlet = bool(sys.argv[3])
        assert isinstance(new_velocity_verlet, bool)
        assert new_velocity_verlet or not new_velocity_verlet, "Safety Measure"
        VELOCITY_VERLET = new_velocity_verlet
    except Exception:
        VELOCITY_VERLET = False
#<Begin: Command Line Args>


#<Begin: MPI Init>
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if not rank:
    statement = '\nUSAGE: mpiexec -n 4 python PATH_TO/mpi_nbody.py (int)TIMESTEPS (int)STATES_TO_SAVE (bool)VELOCITY_VERLET ...'
    if (comm.size == 1) or ('h' in sys.argv) or ('-h' in sys.argv) or ('--h' in sys.argv) or ('help' in sys.argv):
        print statement
        comm.Abort(1)  # exit()
    elif len(sys.argv) == 1:
        print statement
    print '\nUSING: N={:d} TIMESTEPS={:d} STATES_TO_SAVE={:d} VELOCITY_VERLET={:}'.format(comm.size, TIMESTEPS, STATES_TO_SAVE, VELOCITY_VERLET)
#<End: MPI Init>


#<Begin: Physics>
DT = 3.154e8  # Time s
scale_m = 1e5*1.988e30  # Mass kg
scale_x = 3e16  # Length m
scale_v = 400e3/np.sqrt(3.0)  # 1e0*scale_x/DT #Velocity
G = 6.67e-11  # *scale_x*scale_v*scale_v/(scale_m)

# TODO special relativity
# TODO E&M
# TODO collisions or mergers
# TODO spawning
# TODO angular momentum


def random_particle(rank_seed=rank):
    from numpy import pi
    N = 256
    R = 100

    #m = (2*pi-0)*np.random.random_sample((1,))+0
    phi = (2*pi-0)*np.random.random_sample((1,))+0
    costh = (2)*np.random.random_sample((1,))-1
    theta = np.arccos(costh)
    u = np.random.random_sample((1,))

    #<Convert to x,y,z coords>
    r = R * np.power(u, 1/3)
    x = r * np.sin(theta)*np.cos(phi)
    y = r * np.sin(theta)*np.sin(phi)
    z = r * np.cos(theta)
    m = np.random.rand(1)**4.0

    x_y_z_m = np.concatenate([x, y, z, m])
    if np.logical_not(x_y_z_m).sum() > 2:

        return random_particle()

    x_y_z_m[:3] *= scale_x
    x_y_z_m[3] *= scale_m

    vx_vy_vz = ((np.random.rand(3)-0.5)*2.0)**3*scale_v

    ax_ay_az = np.zeros(3)

    return x_y_z_m, vx_vy_vz, ax_ay_az


# THE FOLLOWING LINE MUST BE RUN BEFORE ANY OF THE FOLLOWING CODE
my_xyzm, my_v, my_a = random_particle()


def calculate_cm(many_other_xyzm):
    temp_arr = np.array(many_other_xyzm).T  # TODO shape?
    total_mass = temp_arr[3].sum()
    cmx = np.multiply(temp_arr[0], temp_arr[3]).sum()/total_mass
    cmy = np.multiply(temp_arr[1], temp_arr[3]).sum()/total_mass
    cmz = np.multiply(temp_arr[2], temp_arr[3]).sum()/total_mass
    return cmx, cmy, cmz, total_mass


def calculate_ke(many_other_xyzm, many_other_v):
    masses = np.array(many_other_xyzm).T[3]  # TODO shape?
    # TODO shape??????????????????????????????????????????????????????
    v_squareds = np.square(np.array(many_other_v)).sum(axis=1)
    return np.multiply(masses, v_squareds).sum()/2.0


def calculate_gpe(many_other_xyzm):
    potential = 0.0
    for i in range(1, len(many_other_xyzm)):
        xyzm_i = many_other_xyzm[i]
        potential_i = 0.0
        for xyzm_j in many_other_xyzm[:i]:
            potential_i += xyzm_j[3] / \
                np.sqrt(np.square(xyzm_j[:3]-xyzm_i[:3]).sum())
        potential += xyzm_i[3]*potential_i
    return -1.0*G*potential


def change_in_a_from_xyzm(other_xyzm, xyzm, tolerance=1e-2):
    # TODO
    if (other_xyzm is xyzm) or np.all(other_xyzm == xyzm):
        return np.zeros(3)

    displacement = other_xyzm[:3]-xyzm[:3]

    distance = np.sqrt(np.square(displacement).sum())

    if distance < tolerance:
        return np.zeros(3)  # TODO MERGE OR COLLIDE??

    return G*other_xyzm[3]*displacement/distance**3


def update_state(many_other_xyzm, xyzm=my_xyzm, v=my_v, old_a=my_a):
    a = np.zeros(3)

    if VELOCITY_VERLET:
        # https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
        xyzm[:3] += (v+old_a*DT/2.0)*DT

        # TODO May Want To Do This Destructively
        for other_xyzm in many_other_xyzm:
            a += change_in_a_from_xyzm(other_xyzm, xyzm)

        # a*DT #Use the commented version for Reimann summation; current version is trapezoidal (better)
        v += (a+old_a)*DT/2.0

    else:
        # Trapezoidal Rule, no verlet
        # TODO May Want To Do This Destructively
        for other_xyzm in many_other_xyzm:
            a += change_in_a_from_xyzm(other_xyzm, xyzm)

        # a*DT #Use the commented version for Reimann summation; current version is trapezoidal (better)
        v += (a+old_a)*DT/2.0

        xyzm[:3] += v*DT

    return xyzm, v, a  # TODO May Not Be Needed
#<End: Physics>


#<Begin: Timestep>
if not rank:
    STATES = []
    print '\nINITIAL: my_xyzm', my_xyzm
    rank0_constant = 100/TIMESTEPS
    rank0_window = 99
    rank0_message_base = 'Tstep:{:06d}, Time:{:04.1f}min'
    import time
    rank0_start = time.time()

    def rank0_write(t, many_xyzm=None, many_v=None, many_a=None):
        # Check
        if not TIMESTEPS > t >= 0:
            return
        if not t:
            sys.stdout.write('\n')

        # Wipe
        #sys.stdout.write('\r'+' '*rank0_window)
        # sys.stdout.flush()

        # Generate
        togo = float(TIMESTEPS-t)*(time.time()-rank0_start)/float(t+1)
        rank0_message = rank0_message_base.format(t, togo/60.0)

        if many_xyzm is not None:  # if many_a is None:
            #rank0_message += ', CMx:{:2e}, CMy:{:2e}, CMz:{:2e}, Mtot:{:2e}'
            #CMx,CMy,CMz,Mtot = calculate_cm(many_xyzm)
            # rank0_message = rank0_message.format(CMx,CMy,CMz,Mtot) #calculate_cm(many_xyzm)

            rank0_message += ', GPE:{:2e}'
            gpe = calculate_gpe(many_xyzm)
            rank0_message = rank0_message.format(
                gpe)  # calculate_gpe(many_xyzm)

            if many_v is not None:

                rank0_message += ', KE:{:2e}'
                ke = calculate_ke(many_xyzm, many_v)
                rank0_message = rank0_message.format(
                    ke)  # calculate_ke(many_xyzm, many_v)

                #rank0_message += ', Etot:{:2e}'
                #rank0_message = rank0_message.format(ke+gpe)

        # TODO CM, GPE, KE

        # Write
        if many_xyzm:
            sys.stdout.write('\r'+' '*rank0_window)
            sys.stdout.flush()
        sys.stdout.write('\r'+rank0_message[:rank0_window])
        sys.stdout.flush()

for TIMESTEP in range(TIMESTEPS):
    # TODO Calculate CM for my internal points and do internal nbody

    all_xyzm = comm.allgather(my_xyzm)  # Warning: This Contains Self

    if not TIMESTEP % (TIMESTEPS//STATES_TO_SAVE):
        # TODO: IS THIS IN SOME GUARANTEED ORDER??????????????????????????????????????????
        all_v = comm.gather(my_v, root=0)

        if not rank:
            # TODO Checkpoint
            STATES.append(all_xyzm)  # TODO all_v
            # TODO Readout
            rank0_write(TIMESTEP, all_xyzm, all_v)

    elif not rank:
        rank0_write(TIMESTEP)

    my_xyzm, my_v, my_a = update_state(all_xyzm)

if not rank:
    #sys.stdout.write('\r'+' '*rank0_window+'\n')
    sys.stdout.write('\n')
    sys.stdout.flush()
    print '\nFINAL: my_xyzm', my_xyzm

    import os
    datapath = 'mpi_nbody_out.npy'
    # np.save(os.path.join('SHARED','mpi_nbody_out.npy'),np.array(STATES))
    np.save(datapath, np.array(STATES))
    print '\nData saved as:', datapath, '\n'

#<End: Timestep>
