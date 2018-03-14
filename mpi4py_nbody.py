# ====================================================================
# Author 				: swc21
# Date 					: 2018-03-14 09:41:36
# Project 				: ClusterFiles
# File Name 			: mpi4py_nbody
# Last Modified by 		: swc21
# Last Modified time 	: 2018-03-14 10:55:54
# ====================================================================
#
import numpy as np


def unitize_mass(some_float):
    # mpi4py_nbody
    return SOLAR_MASS*(10**(9*(float(some_float)**3)))
# Particle
# np.random.rand(8)

# First 3 are position, second 3 are velocity, 7 is mass, 8 is unique ID number


def random_particles(n_particles=100):
    # return an Nx8 numpy array of values
    #my_n = 50
    try:
        my_n = abs(int(n_particles))
        assert my_n > 1
    except Exception:
        my_n = 50
    my_particles = np.random.rand(8, my_n)
    assert len(my_particles) == my_n
    assert len(my_particles[0]) == 8
    for my_i, my_particle in enumerate(my_particle):
        my_particle[7] = my_i  # Their last entry is their ID number
