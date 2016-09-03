import logging
import random
import string
import time
#from __future__ import print_function
import h5py as h5
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


import revreaddy as rdy

reload(logging)
logging.basicConfig(
    format='%(asctime)s %(levelname)s: %(message)s',
    datefmt='%Y/%m/%d %I:%M:%S',
    level=logging.DEBUG
)


class Sim_Complex():
    def __init__(self,Diffusion,reactiondistance,micro_reactionrate_forward,micro_reactionrate_backward,micro_reactionrate_complex,alpha,trajectory_length,boxsize,particles):
        """
        order of prameters: first reaction paramters later algorythmic paramter.
        :param boxsize:
        :param Diffusion_S:
        :param Diffusion_P:
        :param reactiondistance:
        :param intrinsikreaction:
        :param length:
        :param particles:
        :param alpha:
        :param Diffusion_E:
        :return:
        """
        self.boxsize = boxsize
        self.Diffusion = Diffusion
        self.reactiondistance = reactiondistance
        self.micro_reactionrate_forward = micro_reactionrate_forward
        self.micro_reactionrate_backward = micro_reactionrate_backward
        self.micro_reactionrate_complex = micro_reactionrate_complex
        self.length = trajectory_length
        self.particles = particles
        self.positions, self.posenzyme = self.envirement()
        self.alpha=alpha

    def stamp(self):
        timestamp = time.strftime("%Y_%m_%d-%H_%M_%S")
        randomstamp = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6))
        stamp = timestamp + "_" + randomstamp
        return stamp

    def envirement(self):
        positions = []
        posenzyme = [0, 0, 0]
        amountenzyme = 0
        while amountenzyme < self.particles:
            # todo periodic boundaryconditions !!!! not anymore couse position of enyzme is 0, 0, 0 .
            pos = np.random.random(3) * (self.boxsize - 0.001) - (self.boxsize * 0.5 - 0.001)
            distance = np.sqrt(
                (posenzyme[0] - pos[0]) ** 2 + (posenzyme[1] - pos[1]) ** 2 + (posenzyme[2] - pos[2]) ** 2)
            if distance > self.reactiondistance:
                amountenzyme += 1
                positions.append(pos)
        positions = np.array(positions)
        return positions, posenzyme

    def observales(self, which_observables):
        observables={}
        for i in range(which_observables):
            if which_observables[i] == "all":
                observables={0 : "MSD",1:"radial",2:"concentration",3:"reactions"}
            elif which_observables == "MSD":
                observables[0] = which_observables[i]
            elif which_observables[i] == "radial":
                observables[1] = which_observables[i]
            elif which_observables[i] == "concentration":
                observables[2] = which_observables[i]
            elif which_observables[i] == "reaction":
                observables[3] = which_observables[i]
        return observables





    def run(self,observables):

        """

        :return: radialdisitbrution fractional
        """
        sim = rdy.Sim("FractionalDiffusion")
        sim.kt = 2.437  # kJ/mol -> ~300 K
        sim.alpha = self.alpha
        sim.timestep = 1.0  # ns
        sim.boxsize = self.boxsize  # nm
        sim.delete_all_particle_types()  # also deletes all interactions etc.
        sim.new_type("S", 0.01, self.Diffusion_S)
        sim.new_type("P", 0.01, self.Diffusion_P)
        sim.new_type("E", 0.01, self.Diffusion_E)
        sim.delete_all_reactions()
        sim.new_enzymatic("E+S<->E+S", 0, 1, 2, self.intrinsikreaction, 0.1, self.reactiondistance)
        sim.delete_all_particles()
        for i in range(self.particles):
            pos = self.positions[i, :]
            sim.add_particle(pos, 0)
            # pos = np.random.random(3) * 9.9 - 4.95
            # sim.add_particle(pos, 1)
            # pos = np.random.random(3) * 9.9 - 4.95
            # sim.add_particle(pos, 2)
        for _ in range(1):
            sim.add_particle(self.posenzyme, 2)
        sim.show_world()
        sim.delete_all_observables()
        mean_e_name = "revdata/" + self.stamp() + "radialdist_e_test1.h5"
        sim.new_radial_distribution(self.length, mean_e_name, np.arange(0., 5., 0.05), np.array([[0, 1]]))
        sim.run(self.length)
        sim.write_observables_to_file()
        sim.delete_all_observables()
        fc = h5.File(mean_e_name, 'r')
        return mean_e_name