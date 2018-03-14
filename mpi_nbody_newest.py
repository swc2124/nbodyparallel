# ====================================================================
# Author 				: swc21
# Date 					: 2018-03-14 09:42:27
# Project 				: ClusterFiles
# File Name 			: mpi_nbody_newest
# Last Modified by 		: swc21
# Last Modified time 	: 2018-03-14 10:29:53
# ====================================================================
#
#<Begin: Imports>
from mpi4py import MPI
import numpy as np
import sys
import operator
import os
import time
import mpi_nbody_newest_lib as library  # Homemade Code!
#<End: Imports>
#<Begin: Command Line Args>
BIG_N = 100000
TIMESTEPS = 100
if len(sys.argv) > 1:
    try:
        new_timesteps = int(sys.argv[1])
        assert isinstance(new_timesteps, int)
        assert BIG_N > new_timesteps >= 2, "Safety Measure"
        TIMESTEPS = new_timesteps
    except Exception:
        TIMESTEPS = 100
STATES_TO_SAVE = 10
if len(sys.argv) > 2:
    try:
        new_states_to_save = int(sys.argv[2])
        assert isinstance(new_states_to_save, int)
        assert min(BIG_N, TIMESTEPS) >= new_states_to_save >= 2, "Safety Measure"
        STATES_TO_SAVE = new_states_to_save
    except Exception:
        STATES_TO_SAVE = 10
'''
VELOCITY_VERLET = False
if len(sys.argv) > 3:
    try:
        new_velocity_verlet = bool(sys.argv[3])
        assert isinstance(new_velocity_verlet, bool)
        assert new_velocity_verlet or not new_velocity_verlet, "Safety Measure"
        VELOCITY_VERLET = new_velocity_verlet
    except Exception:
        VELOCITY_VERLET = False
'''
#<Begin: Command Line Args>
#<Begin: MPI Init>
comm_world = MPI.COMM_WORLD
group_world = comm_world.group
rank_world = comm_world.Get_rank()
if not rank_world:
    def partition(list_of_min_max_tuples, n_employees=comm_world.size):
        for listing in list_of_min_max_tuples:
            assert listing[0] >= 0
            assert listing[1] is None or listing[1] >= listing[0]
        departments = []
        while len(departments) < len(list_of_min_max_tuples):
            departments.append([])
        employees = range(n_employees)
        employees.reverse()
        for department, listing in zip(departments, list_of_min_max_tuples):
            while len(department) < listing[0]:  # Min employees per department
                # This is allowed to fail: means you have too few n
                department.append(employees.pop())
        for department, listing in zip(departments, list_of_min_max_tuples):
            # Max employees per department
            while listing[1] is None or len(department) < listing[1]:
                try:
                    # This is allowed to fail: means you have too few n
                    department.append(employees.pop())
                except Exception:
                    break
        if not rank_world and len(employees):
            print 'Unused processes that will not get tasked:', employees
        return departments
    min_observers, max_observers = 1, 1
    min_workers, max_workers = 2, None  # Will take any additional ranks as workers
    # min_others,max_others = 0,0 #Placeholder
    observers, workers = partition(
        [(min_observers, max_observers), (min_workers, max_workers)])
    managers = [observers[0], workers[0]]
else:
    observers, workers, managers = None, None, None
observers = comm_world.bcast(observers)
workers = comm_world.bcast(workers)
managers = comm_world.bcast(managers)
#head_observer = observers[0]
#head_worker = workers[0]
#head_other = others[0]


def am_manager(r=rank_world):
    return bool(r in managers)


def am_observer(r=rank_world):
    return bool(r in observers)


def am_head_observer(r=rank_world):
    return am_observer(r) and am_manager(r)  # bool(r == head_observer)


def am_worker(r=rank_world):
    return bool(r in workers)


def am_head_worker(r=rank_world):
    return am_worker(r) and am_manager(r)  # bool(r == head_worker)


if am_head_observer():
    # print comm_world.gather(rank_world)
    # (bool)VELOCITY_VERLET ...'
    statement = '\nUSAGE: mpiexec -n 4 python PATH_TO/mpi_nbody.py (int)TIMESTEPS (int)STATES_TO_SAVE'
    if ('h' in sys.argv) or ('-h' in sys.argv) or ('--h' in sys.argv) or ('help' in sys.argv):
        print statement
        comm_world.Abort(1)  # exit()
    elif len(sys.argv) == 1:
        print statement
    # VELOCITY_VERLET={:}'
    varstring = '\nUSING: Nbodies={:d} TIMESTEPS={:d} STATES_TO_SAVE={:d}'
    # , VELOCITY_VERLET)
    print varstring.format(len(workers), TIMESTEPS, STATES_TO_SAVE)
group_observers = group_world.Incl(observers)
comm_observers = comm_world.Create(group=group_observers)
group_workers = group_world.Incl(workers)
comm_workers = comm_world.Create(group=group_workers)
group_managers = group_world.Incl(managers)
comm_managers = comm_world.Create(group=group_managers)
'''
if am_observer():
    rank_observers = comm_observers.Get_rank()
elif am_worker():
    rank_workers = comm_workers.Get_rank()
if am_manager():
    rank_managers = comm_managers.Get_rank()
'''
comm_world.barrier()
#<End: MPI Init>
try:
    #<Begin: Main>
    if am_head_observer():
        rank0_window = 70
        rank0_message_base = 'Tstep:{:06d}, Time:{:04.1f}min'
        rank0_start = time.time()
        my_States = library.States()
        KE_0, GPE_0 = None, None

        def observe(t, state=None):
            # Check
            if not TIMESTEPS > t >= 0:
                return
            if not t:
                sys.stdout.write('\n')
            global KE_0, GPE_0
            # Generate
            togo = float(TIMESTEPS-t)*(time.time()-rank0_start)/float(t+1)
            if t+1 == TIMESTEPS:
                togo = time.time()-rank0_start
            rank0_message = rank0_message_base.format(t, togo/60.0)
            if state is not None:
                assert isinstance(state, library.State)
                sys.stdout.write('\r'+' '*rank0_window)
                sys.stdout.flush()
                # TODO-------------------------------------------------------------------------------
                rank0_message += ', KE:{:2e}'
                ke = state.KE()
                if KE_0 is None:
                    KE_0 = ke
                    rank0_message = rank0_message.format(ke)
                else:
                    rank0_message = rank0_message.format(ke/KE_0)
                # TODO-------------------------------------------------------------------------------
                rank0_message += ', GPE:{:2e}'
                gpe = state.GPE()
                if GPE_0 is None:
                    GPE_0 = gpe
                    rank0_message = rank0_message.format(gpe)
                else:
                    rank0_message = rank0_message.format(gpe/GPE_0)
                # TODO-------------------------------------------------------------------------------
                rank0_message += ', SUM:{:2e}'
                value = gpe+2*ke
                rank0_message = rank0_message.format(value)
                # TODO-------------------------------------------------------------------------------
                my_States.add(state)
            sys.stdout.write('\r'+rank0_message[:rank0_window])
            sys.stdout.flush()
    elif am_worker():
        # THE FOLLOWING LINE MUST BE RUN BEFORE ANY OF THE FOLLOWING CODE
        #my_xyzm, my_v, my_a = random_particle()
        '''
        if am_head_worker():
                all_xyzm, all_v = library.add_particles(comm_workers.size)  
        else:
                all_xyzm, all_v = None,None
        my_xyzm = comm_workers.scatter(sendobj=all_xyzm, root=head_worker) #recvobj
        my_xyzm = np.asarray(my_xyzm)
        my_v = comm_workers.scatter(sendobj=all_v, root=head_worker)
        my_v = np.asarray(my_v)
        my_a = np.zeros(3)
        #del all_xyzm
        #del all_v
        '''
        comm_workers.barrier()
        if am_head_worker():
            # Oldest
            #all_Particles = map(library.Particle, library.add_particles(comm_workers.size))
            # Old
            #junk = library.add_particles(comm_workers.size)
            #all_Particles = [library.Particle(a,b) for a,b in zip(junk[0], junk[1])]
            # New
            all_Particles = library.add_particles(comm_workers.size)
        else:
            all_Particles = None
        comm_workers.barrier()
        #my_Particle = None
        my_Particle = comm_workers.scatter(all_Particles)
        #my_Particle = comm_workers.Scatter(all_Particles, [library.Particle], root=head_worker)
    for TIMESTEP in range(TIMESTEPS):
        #<Begin: Timestep>
        comm_world.barrier()
        # print 't:'+str(TIMESTEP)+' r_w:'+str(rank_world)
        if am_worker():
            my_xyzm = my_Particle.get_xyzm()
            all_xyzm = comm_workers.allgather(my_xyzm)
            new_xyzm, new_v, new_a = library.update_state(
                all_xyzm, my_xyzm, my_Particle.vx_vy_vz.copy(), my_Particle.ax_ay_az.copy())
            my_Particle = library.Particle(new_xyzm, new_v, new_a)
        if not TIMESTEP % (TIMESTEPS//STATES_TO_SAVE):
            if am_worker():
                STATE = library.State(
                    TIMESTEP, comm_workers.gather(my_Particle))
                if am_head_worker():
                    #comm_world.send(STATE, dest=0, tag=555)
                    comm_managers.send(STATE, dest=int(
                        not comm_managers.rank), tag=555)
            elif am_head_observer():
                #STATE = comm_world.recv(source=1, tag=555)
                STATE = comm_managers.recv(
                    source=int(not comm_managers.rank), tag=555)
                observe(TIMESTEP, STATE)
        elif am_head_observer():
            observe(TIMESTEP)
        # print 't:'+str(TIMESTEP)+' r_w:'+str(rank_world)
        comm_world.barrier()
        #<End: Timestep>
    if am_head_observer():
        sys.stdout.write('\n')
        sys.stdout.flush()
        datapath = 'mpi_nbody_out.npy'
        my_States.save(datapath)
        print '\nData saved as:', datapath
    # TODO More in Main ?? Analysis ?? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if am_head_observer():
        print '\nDone!\n'
    #<End: Main>
except Exception as err:
    try:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        print rank_world, exc_tb.tb_lineno, exc_type, err
    except Exception:
        print 'Bad traceback! No Further info available!'
    finally:
        comm_world.Abort(1)
finally:
    # my_States.close() ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pass
