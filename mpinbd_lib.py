# ====================================================================
# Author 				: swc21
# Date 					: 2018-03-14 09:41:11
# Project 				: ClusterFiles
# File Name 			: mpinbd_lib
# Last Modified by 		: swc21
# Last Modified time 	: 2018-03-14 10:29:04
# ====================================================================
#
# version 2 mpi nbody lib
#<Begin: Parameters>
TIMESTEPS = 5
STATES_TO_SAVE = 1
#<End: Parameters>
#<Begin: Imports>
import numpy as np
import operator
import os
import sys
import time

from mpi4py import MPI
#<End: Imports>
#<Begin: Constants>
DT = 3e-2*3.154e7  # Time year in s
C = 2.99*1e8  # m/s
C_squared = C**2  # (m/s)**2
G = 6.67e-11  # N*(m**2)/kg**2_
m_e = 9.109e-31  # kg_mass of electron
m_p = 1.6e-27  # kg_protonmass
u_0 = 4.0e-7*np.pi  # N/amp**2_permeability
E_0 = 1.0/(C_squared*u_0)  # Farad/m_permitivity
K_e = 1.0/4*np.pi*E_0  # N*(m**2)/Coulombs**2_Coulombsconstant
scale_m = 1e-1*1.988e30  # Mass kg
scale_x = 3e16  # Length m
scale_v = 4e2  # m/s_Velocity
scale_q = 1.6e-19  # protoncharge_Coulombs#look in lib
#<End: Constants>
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
            # TODO:kill extra workers
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
#<Begin: Class Definitions> #TODO Will be different for every simulation!
dScalar = np.float64
dVector = np.dtype([('x', dScalar), ('y', dScalar), ('z', dScalar)])
dTensor = np.dtype([('x', dVector), ('y', dVector), ('z', dVector)])
dParticle = np.dtype([('time', dScalar), ('mass', dScalar), ('charge', dScalar), ('radius', dScalar), ('position', dVector),
                      ('velocity', dVector), ('acceleration', dVector), ('angular_momentum', dVector), ('moment_of_Inertia', dTensor)])
#dState_Expanded = np.dtype(('time', dScalar), ('particles', dParticle, comm_workers.size))
#dState_ = np.dtype(('time', dScalar), ('particle', dParticle))
# TODO class Particle


def Acceleration(other, one):
    assert other.dtype is dParticle
    assert one.dtype is dParticle
    assert other.shape == (comm_workers.size,)
    assert one.shape == (1,)
    # TODO
    return np.zeros((1,), dtype=dVector)


class State:
    def __init__(self, one_Particle_or_State):
        if type(one_Particle_or_State) is State:
            self.myTime = one_Particle_or_State.myTime
            self.myParticle = one_Particle_or_State.myParticle
            self.myParticles = one_Particle_or_State.myParticles
            self.isExpanded = one_Particle_or_State.isExpanded
        else:
            assert one_Particle_or_State.dtype == dParticle
            #assert one_Particle_or_State.shape == (1,)
            # TODO More assertions
            self.myTime = 0.0
            self.myParticle = one_Particle_or_State.copy().view(np.recarray)
            self.myParticles = None
            self.isExpanded = False

    def step_expand(self):
        # comm_workers.barrier()
        print 'pre-'+str(lib.rank_world)
        lib.comm_workers.barrier()
        assert not self.isExpanded
        self.myParticles = comm_workers.allgather(self.myParticle)
        self.myParticle = None
        self.isExpanded = True
        print 'post-'+str(lib.rank_world)

    def step_compress(self, dt=DT):
        # comm_workers.barrier()
        assert self.isExpanded
        # TODO, remove my particle from self.Particles?
        self.__update__(dt)  # This bring myParticle up-to-date
        self.myParticles = None
        self.isExpanded = False
        self.myTime += dt

    def __update__(self, dt):
        # This bring myParticle up-to-date
        myParticle_old = myParticle.copy()
        myParticle.acceleration = np.zeros((1,), dtype=dVector)
        for otherParticle in myParticles:
            if myParticle != otherParticle:
                myParticle.acceleration += Acceleration(
                    myParticle_old.copy(), otherParticle.copy())
        # TODO MORE COMPLEX INTEGRATORS
        myParticle.velocity += (myParticle.acceleration +
                                myParticle_old.acceleration)*dt/2.0
        myParticle.position += (myParticle.velocity +
                                myParticle_old.velocity)*dt/2.0
#<End: Class Definitions>
#<Begin: Factory Functions>


def add_particles(N):
    # TODO
    return np.zeros((N), dtype=dParticle)
#<End: Factory Functions>
