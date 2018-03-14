# ====================================================================
# Author 				: swc21
# Date 					: 2018-03-14 09:42:27
# Project 				: ClusterFiles
# File Name 			: mpi_nbody_newest_lib
# Last Modified by 		: swc21
# Last Modified time 	: 2018-03-14 10:28:54
# ====================================================================
#
#<Begin: Imports>
import numpy as np
import operator
import random
import math
import Solersystem2 as sol
#<End: Imports>
#<Begin: Physics>
'''
N_Stinkers = 1
DT = 1000*3.154e7 # Time year in s
C = 2.99*1e8 # m/s
C_squared = C**2 # (m/s)**2
sol.G = 6.67e-11 #N*(m**2)/kg**2_ 
m_e = 9.109e-31 #kg_mass of electron
m_p = 1.6e-27 # kg_protonmass
u_0 = 4.0e-7*np.pi #N/amp**2_permeability
E_0 = 1.0/(C_squared*u_0) #  Farad/m_permitivity 
K_e = 1.0/4*np.pi*E_0 # N*(m**2)/Coulombs**2_Coulombsconstant
sol.scale_m = 1e-1*1.988e30 # Mass kg
sol.scale_x = 3e16 # Length m 
scale_v = 4e2 # m/s_Velocity
sol.scale_q = 1.6e-19 #protoncharge_Coulombs
'''
# from Solersystem2 import N_Stinkers, C, C_squared, DT, sol.scale_m, sol.scale_x, scale_v, sol.scale_q, sol.G
#from Solersystem2 import add_particles as r_v_a_m
# sol.scale_q = 1.602e-19 #Coulombs
SPECIAL_RELATIVISTIC = True
# VELOCITY_VERLET = False  #(Don't set this to True without good reason!)
# TODO E&M
# TODO collisions or mergers
# TODO spawning
# TODO angular momentum


class Particle:
    def __init__(self, x_y_z_m_q, vx_vy_vz=np.zeros(3), ax_ay_az=np.zeros(3), mass=None, charge=None, radius=0.0, Lx_Ly_Lz=np.zeros(3), bx_by_bz=np.zeros(3)):
        # DO NOT CHANGE THE ABOVE ORDER #FOR BACKWARDS COMPATABILITY
        x_y_z_m_q = np.array(x_y_z_m_q).astype(float)
        if x_y_z_m_q.shape == (3,):
            self.x_y_z = x_y_z_m_q[:3]
            if mass is not None:
                self.m = float(mass)
            else:
                self.m = float(sol.scale_m)
            if charge is not None:
                self.q = float(charge)
            else:
                self.q = float(sol.scale_q)
        elif x_y_z_m_q.shape == (4,):
            self.x_y_z = x_y_z_m_q[:3]
            self.m = float(x_y_z_m_q[3])
            if charge is not None:
                self.q = float(charge)
            else:
                self.q = float(sol.scale_q)
        elif x_y_z_m_q.shape == (5,):
            self.x_y_z = x_y_z_m_q[:3]
            self.m = float(x_y_z_m_q[3])
            self.q = float(x_y_z_m_q[4])
        else:
            raise AssertionError('x_y_z_m_q: '+str(x_y_z_m_q))
        assert self.m >= 0.0
        vx_vy_vz = np.array(vx_vy_vz).astype(float)
        assert vx_vy_vz.shape == (3,)
        self.vx_vy_vz = vx_vy_vz
        ax_ay_az = np.array(ax_ay_az).astype(float)
        assert ax_ay_az.shape == (3,)
        self.ax_ay_az = ax_ay_az
        assert Lx_Ly_Lz.shape == (3,)
        self.Lx_Ly_Lz = Lx_Ly_Lz.astype(float)
        assert bx_by_bz.shape == (3,)
        self.bx_by_bz = bx_by_bz.astype(float)
        assert isinstance(float(radius), float)
        self.r = float(radius)
        assert self.r >= 0.0
        self.is_Massless = not self.m
        self.is_Chargeless = not self.q
        self.is_Point = not self.r
        # TODO More

    def get_xyzm(self):
        temp = np.zeros(4, dtype=float)
        temp[:3] = self.x_y_z
        temp[3] = self.m
        return temp

    def set_xyzm(self, temp):
        assert temp.shape == (4,)
        self.x_y_z = temp[:3].astype(float)
        self.m = float(temp[3])

    def get_xyzq(self):
        temp = np.zeros(4, dtype=float)
        temp[:3] = self.x_y_z
        temp[3] = self.q
        return temp

    def set_xyzq(self, temp):
        assert temp.shape == (4,)
        self.x_y_z = temp[:3].astype(float)
        self.q = float(temp[3])

    def get_xyzmq(self):
        temp = np.zeros(5, dtype=float)
        temp[:3] = self.x_y_z
        temp[3] = self.m
        temp[4] = self.q
        return temp

    def set_xyzmq(self, temp):
        assert temp.shape == (5,)
        self.x_y_z = temp[:3].astype(float)
        self.m = float(temp[3])
        self.q = float(temp[4])


class Field:
    def __init__(self):
        # TODO
        pass


class State:
    def __init__(self, time, zero_or_more_Particle_objects=None, zero_or_more_Field_objects=None):
        self.Time = float(time)
        self.Particles = []
        if zero_or_more_Particle_objects is not None:
            for p in zero_or_more_Particle_objects:
                assert isinstance(p, Particle)
                self.Particles.append(p)
        self.Fields = []
        if zero_or_more_Field_objects is not None:
            for f in zero_or_more_Field_objects:
                assert isinstance(f, Field)
                self.Fields.append(f)

    def GPE(self):
        return calculate_gpe(map(Particle.get_xyzm, self.Particles))

    def KE(self):
        return calculate_ke(map(Particle.get_xyzm, self.Particles), map(operator.attrgetter('vx_vy_vz'), self.Particles))


class States:
    def __init__(self, zero_or_more_State_objects=None, cache=False):
        self.States = []
        self.is_Cached = bool(cache)
        if zero_or_more_State_objects is not None:
            for s in zero_or_more_State_objects:
                assert isinstance(s, State)
                # self.States.append(s)
                self.add(s)
        # self.States.sort(key=operator.attrgetter('Time'))

    def __iter__(self):
        return self.States

    def extend(self, States_object=None):
        if isinstance(States_object, States):
            # TODO Destructive??
            for state in States_object:
                self.add(state)
        else:
            for s in States_object:
                assert isinstance(s, State)
                self.add(s)

    def add(self, one_State_object=None):
        assert isinstance(one_State_object, State)
        if self.is_Cached:
            pass
        else:
            self.States.append(one_State_object)

    def save(self, path='mpi_nbody_newest_out'):
        if not self.is_Cached:
            np.save(path, np.asarray(
                [[p.get_xyzm() for p in s.Particles] for s in self.States]))
        else:
            pass

    def load(self, path='mpi_nbody_newest_out'):
        print 'load not avilable yet!'
        # TODO
        return
        if not self.is_Cached:
            self.extend(np.load(path))
        else:
            pass


def add_particles(N):
    return [Particle(*sol_thing) for sol_thing in sol.add_particles(N)]
    #[0., 0., 0.]


def change_in_a_from_xyzm(other_xyzm, xyzm, tolerance=1e6):
    if (other_xyzm is xyzm) or np.all(other_xyzm == xyzm):
        return np.zeros(3)
    displacement = other_xyzm[:3]-xyzm[:3]
    distance = np.sqrt(np.square(displacement).sum())
    if distance < tolerance:
        return np.zeros(3)  # TODO MERGE OR COLLIDE??
    return sol.G*other_xyzm[3]*displacement/distance**3


def gamma(velocity):
    assert isinstance(velocity, np.ndarray)
    assert velocity.shape == (3,)
    return float(1.0/np.sqrt(1.0+np.square(velocity).sum()/sol.C_squared))


def sr_helper2(_v, _i=np.zeros(3)):
    assert isinstance(_v, np.ndarray)
    assert _v.shape == (3,)
    assert isinstance(_i, np.ndarray)
    assert _i.shape == (3,)
    return (gamma(_v)*_v)+_i


def new_v_from_impulse_and_old_v(old_v, impulse):
    # impulse = (a_effective)*dt #Must have units of velocity
    # a_effective = (new_a+old_a)/2.0 #This is a 2nd Order Method
    return sr_helper2(sr_helper2(old_v, impulse))


def update_state(many_other_xyzm, xyzm, v, old_a):
    a = np.zeros(3)
    # if VELOCITY_VERLET:
    # https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
    #xyzm[:3] += (v+old_a*DT/2.0)*DT
    # TODO May Want To Do This Destructively
    for other_xyzm in many_other_xyzm:
        a += change_in_a_from_xyzm(other_xyzm, xyzm)
    impulse = xyzm[3]*(a+old_a)*sol.DT/2.0  # a*DT
    if SPECIAL_RELATIVISTIC:
        v = new_v_from_impulse_and_old_v(v.copy(), impulse.copy())
    else:
        v += impulse
    # if not VELOCITY_VERLET:
    # 2nd order integrator
    xyzm[:3] += v*sol.DT
    # print 'us:'+str(rank_workers)
    return xyzm, v, a  # TODO May Not Be Needed
# Physics for all


def calculate_cm(many_other_xyzm):
    temp_arr = np.array(many_other_xyzm).T  # TODO shape?
    total_mass = temp_arr[3].sum()
    cmx = np.multiply(temp_arr[0], temp_arr[3]).sum()/total_mass
    cmy = np.multiply(temp_arr[1], temp_arr[3]).sum()/total_mass
    cmz = np.multiply(temp_arr[2], temp_arr[3]).sum()/total_mass
    return cmx, cmy, cmz, total_mass


def calculate_ke(many_other_xyzm, many_other_v):
    masses = np.array(many_other_xyzm).T[3]  # TODO shape?
    v_squareds = np.square(np.array(many_other_v)).sum(axis=1)  # TODO shape?
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
    return -1.0*sol.G*potential
#<End: Physics>
