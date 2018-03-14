# ====================================================================
# Author 				: swc21
# Date 					: 2018-03-14 09:41:11
# Project 				: ClusterFiles
# File Name 			: mpinbd
# Last Modified by 		: swc21
# Last Modified time 	: 2018-03-14 10:29:02
# ====================================================================
#
# Version 2 mpi nbody
import sys
import mpinbd_lib as lib  # Homemade Code!
#<Begin: loading>
if lib.am_head_observer():
    print "Loading..."
    '''
	#print comm_world.gather(rank_world)
	statement = '\nUSAGE: mpiexec -n 4 python PATH_TO/mpi_nbody.py (int)TIMESTEPS (int)STATES_TO_SAVE'# (bool)VELOCITY_VERLET ...'
	if ('h' in sys.argv) or ('-h' in sys.argv) or ('--h' in sys.argv) or ('help' in sys.argv):
		print statement
		comm_world.Abort(1) #exit()
	elif len(sys.argv)==1:
		print statement
	varstring='\nUSING: Nbodies={:d} TIMESTEPS={:d} STATES_TO_SAVE={:d}'# VELOCITY_VERLET={:}'
	print varstring.format(len(workers), TIMESTEPS, STATES_TO_SAVE)#, VELOCITY_VERLET)
	'''
#<End: loading>
#<Begin: Main>
try:
    #<step 1>
    if lib.am_head_observer():
        sys.stdout.write('\n')
        import time
        rank0_window = 70
        rank0_message_base = 'Tstep:{:06d}, Time:{:04.1f}min'
        rank0_start = time.time()

        def observe(t, state=None):
            # Check
            print "observe()", lib.rank_world
            if not isinstance(t, int) or not lib.TIMESTEPS > t >= 0:
                return
            if t+1 == lib.TIMESTEPS:
                togo = time.time()-rank0_start
            else:
                togo = float(lib.TIMESTEPS-t) * \
                    (time.time()-rank0_start)/float(t+1)
            rank0_message = rank0_message_base.format(t, togo/60.0)
            # TODO:
            if state is not None and type(state) is lib.State:
                pass
            sys.stdout.write('\r'+rank0_message[:rank0_window])
            sys.stdout.flush()
    elif lib.am_worker():
        # The goal of these steps is to have a single particle to contribute to the global state
        lib.comm_workers.barrier()
        if lib.am_head_worker():
            all_Particles = lib.add_particles(lib.comm_workers.size)
        else:
            all_Particles = None
        lib.comm_workers.barrier()
        my_Particle = lib.comm_workers.scatter(all_Particles)
        lib.comm_workers.barrier()
        # In the previous steps, each worker just somehow needs to make a Particle
        my_State = lib.State(my_Particle)
        lib.comm_workers.barrier()
    #<step 2>
    for TIMESTEP in range(lib.TIMESTEPS):
        # WARNING, THE ORDER OF THESE IF BLOCKS REALLY MATTERS A LOT! DO NOT MOVE!
        if lib.am_worker():
            print 'pre-my_State.step_expand()'+str(lib.rank_world)
            lib.comm_workers.barrier()
            State.step_expand(my_State)
            # my_State.step_expand() #lib.comm_workers.barrier()
            print 'post-my_State.step_expand()'+str(lib.rank_world)
        if not TIMESTEP % (lib.TIMESTEPS//lib.STATES_TO_SAVE):
            if lib.am_head_worker():
                # TODO Cannot have more than one other manager
                lib.comm_managers.send(lib.State(my_State), dest=int(
                    not lib.comm_managers.rank), tag=TIMESTEP)  # tag=555)
            elif lib.am_head_observer():
                # TODO Cannot have more than one other manager
                observe(TIMESTEP, lib.comm_managers.recv(source=int(
                    not lib.comm_managers.rank), tag=TIMESTEP))  # tag=555))
        elif am_head_observer():
            observe(TIMESTEP)
        if lib.am_worker():
            print 'pre-my_State.step_compress()'+str(lib.rank_world)
            lib.comm_workers.barrier()
            my_State.step_compress()  # dt=lib.DT #lib.comm_workers.barrier()
            print 'post-my_State.step_compress()'+str(lib.rank_world)
    #<step 3>
    if am_head_observer():
        print "am_head_observer()", lib.rank_world
        sys.stdout.write('\n')
        sys.stdout.flush()
        #datapath = 'mpi_nbody_out.npy'
        # my_States.save(datapath)
        # print '\nData saved as:',datapath
        print '\nDone!\n'
except Exception as err:
    try:
        exc_type, exc_obj, exc_tb = sys.exc_info()
        print rank_world, exc_tb.tb_lineno, exc_type, err
    except Exception:
        print 'Bad traceback! No Further info available!'
    finally:
        lib.comm_world.Abort(1)
finally:
    # TODO open_files.close() ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pass
#<End: Main>
