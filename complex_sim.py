import logging
import os
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
    def __init__(self,Diffusion,reactiondistance,micro_reactionrate_forward,micro_reactionrate_backward,micro_reactionrate_complex,alpha,trajectory_length,boxsize,particles,timestep):
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
        self.observables={}
        self.timestep=timestep
        self.path=self.path()

    def info(self):
        name='Info'+'.txt'
        with open(self.path+name, 'w') as info:
            info.write("This is the description of the setup of the Simulation")
            info.write("boxsize= %d" %self.boxsize)
            info.write("Diffusion= %d" %self.Diffusion)
            info.write("reactiondistance= %d" %self.reactiondistance)
            info.write("reactionrate_f= %d" %self.micro_reactionrate_forward)
            info.write("reactionrate_b= %d" %self.micro_reactionrate_backward)
            info.write("reactionrate_c= %d" %self.micro_reactionrate_complex)
            info.write("alpha= %d" %self.alpha)
            info.write("trajectory length= %d" %self.length)
            info.write("amount of particles= %d" %self.particles)
            info.write("timestep=%d" %self.timestep)

    def stamp(self):
        timestamp = time.strftime("%Y_%m_%d-%H_%M_%S")
        randomstamp = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(6))
        stamp = timestamp + "_" + randomstamp
        return stamp
    def path(self):
        path= "/srv/data/Janekg89/revdata/gute/" +"rate_f"+str(self.micro_reactionrate_forward)+"rate_b"+str(self.micro_reactionrate_backward)+"rate_c"+str(self.micro_reactionrate_complex)+"alpha"+str(self.alpha)+"S_num"+str(self.particles)+"length"+str(self.length)+"R"+str(self.reactiondistance)+"boxsize"+str(self.boxsize)+"D"+str(self.Diffusion)+"/"
        if not os.path.exists(path):
            os.makedirs(path)
        return path

    def envirement(self):
        positions = []
        posenzyme = [0, 0, 0]
        amountenzyme = 0
        while amountenzyme < self.particles:
            # todo periodic boundaryconditions !!!! not anymore couse position of enyzme is 0, 0, 0 .
            pos = np.random.random(3) * (self.boxsize - 0.00001) - (self.boxsize * 0.5 - 0.00001)
            distance = np.sqrt(
                (posenzyme[0] - pos[0]) ** 2 + (posenzyme[1] - pos[1]) ** 2 + (posenzyme[2] - pos[2]) ** 2)
            #if distance > self.reactiondistance:
            amountenzyme += 1
            positions.append(pos)
        positions = np.array(positions)
        return positions, posenzyme

    def observales(self, which_observables):

        for i in range(which_observables):
            if which_observables[i] == "all":
                self.observables={0 : "MSD",1:"radial",2:"concentration",3:"number_of_particle"}
            elif which_observables == "MSD":
                self.observables[0] = which_observables[i]
            elif which_observables[i] == "radial":
                self.observables[1] = which_observables[i]
            elif which_observables[i] == "number_of_particle":
                self.observables[2] = which_observables[i]
            elif which_observables[i] == "reaction":
                self.observables[3] = which_observables[i]


    def run(self):

        """

        runs the simulation
        """

        if self.alpha==1.0:
            sim = rdy.Sim("")
        else:
            sim=rdy.Sim("FractionalDiffusion")

        sim.kt = 1  # kJ/mol -> ~300 K
        sim.alpha = self.alpha
        sim.timestep = self.timestep  # ns
        sim.boxsize = self.boxsize  # nm
        sim.delete_all_particle_types()  # also deletes all interactions etc.
        sim.new_type("S", 0.001, self.Diffusion)
        sim.new_type("P", 0.001, self.Diffusion)
        sim.new_type("E", 0.001, 0)
        sim.new_type("C",0.001,0)
        sim.delete_all_reactions()
        sim.new_fusion("S+E<->C",0,2,3,self.micro_reactionrate_forward,self.micro_reactionrate_backward,self.reactiondistance)
        sim.configure_fusion(0,[],1.,0.5*self.reactiondistance,1,1,1,0,1)
        sim.new_fusion("P+E<-C", 1,2,3,0, self.micro_reactionrate_complex, self.reactiondistance)
        sim.configure_fusion(1,[],1.,0.5*self.reactiondistance,1,1,1,0,1)
        sim.delete_all_particles()
        for i in range(self.particles):
            pos = self.positions[i, :]
            sim.add_particle(pos, 0)
        for _ in range(1):
            sim.add_particle(self.posenzyme, 2)
        #sim.show_world()
        sim.delete_all_observables()

        for ii in self.observables:
            if ii == "all":
                combistamp=self.stamp()

                msd=self.path + combistamp + "msd.h5"
                sim.new_mean_squared_displacement(1,msd,0)

                S_and_E_to_C=self.path + combistamp + "reac_num_S_and_E_to_C.h5"
                P_and_E_to_C=self.path + combistamp + "reac_num_S_and_P_to_C.h5"

                sim.new_reaction_counter(1,S_and_E_to_C,"S+E<->C")
                sim.new_reaction_counter(1,P_and_E_to_C,"P+E<-C")

                radial_S = self.path + combistamp + "radialS.h5"
                radial_P = self.path + combistamp + "radialP.h5"
                sim.new_radial_distribution(self.length, radial_S, np.arange(0, (self.boxsize*1.0)/2, ((self.boxsize*1.0)/(2.0*30.0) )), np.array([[0,2],[0,3]]))
                sim.new_radial_distribution(self.length, radial_P, np.arange(0, (self.boxsize*1.0)/2, ((self.boxsize*1.0)/(2.0*30.0) )), np.array([[1, 2],[1,3]]))


                number_S=self.path + combistamp + "number_s.h5"
                number_P=self.path + combistamp + "number_p.h5"
                number_C=self.path + combistamp + "number_c.h5"

                sim.new_particle_numbers(1,number_S,0)
                sim.new_particle_numbers(1,number_P,1)
                sim.new_particle_numbers(1,number_C,3)

            elif  ii == "MSD":
                msd=self.path + self.stamp() + "msd.h5"
                sim.new_mean_squared_displacement(1,msd,0)

            elif ii == "radial":
                """
                gives the radial distribution of S around E and C and of P around E and C on the last frame of the simulation

                r_momentan=1
                r_radial=[r_momentan]
                v=self.boxsize**3/50.0
                for _ in range(30):

                    dr=-r_momentan+(r_momentan**3+((v*3)/(4*np.pi)))**(1.0/3.0)
                    r_momentan=r_momentan+dr
                    r_radial.append(r_momentan)
                """
                combistamp=self.stamp()
                radial_S = self.path + combistamp + "radialS.h5"
                radial_P = self.path + combistamp + "radialP.h5"
                number_S=self.path + combistamp + "number_s.h5"
                number_P=self.path + combistamp + "number_p.h5"
                number_C=self.path + combistamp + "number_c.h5"
                sim.new_particle_numbers(1,number_S,0)
                sim.new_particle_numbers(1,number_P,1)
                sim.new_particle_numbers(1,number_C,3)
                sim.new_radial_distribution(self.length, radial_S, np.arange(0, self.boxsize/2, (self.boxsize/(2.0*30.0) )), np.array([[0,2],[0,3]]))
                sim.new_radial_distribution(self.length, radial_P, np.arange(0, self.boxsize/2, (self.boxsize/(2.0*30.0) )), np.array([[1, 2],[1,3]]))
                #sim.new_radial_distribution(self.length, radial_S,r_radial, np.array([[2,0],[3,0]]))
                #sim.new_radial_distribution(self.length, radial_P,r_radial, np.array([[2, 1],[3,1]]))

            elif ii == "number_of_particle":
                """
                gives number of particles  of s,p and c over time
                """
                for iii in self.observables:
                    if iii == "radial":
                        pass
                    else:
                        number_S=self.path + self.stamp() + "number_s.h5"
                        number_P=self.path + self.stamp() + "number_p.h5"
                        number_C=self.path + self.stamp() + "number_c.h5"
                        sim.new_particle_numbers(1,number_S,0)
                        sim.new_particle_numbers(1,number_P,1)
                        sim.new_particle_numbers(1,number_C,3)
            elif ii == "reaction":
                """
                gives number of reactions over time
                """
                S_and_E_to_C=self.path + self.stamp() + "reac_num_S_and_E_to_C.h5"
                P_and_E_to_C=self.path + self.stamp() + "reac_num_S_and_P_to_C.h5"
                sim.new_reaction_counter(1,S_and_E_to_C,"S+E<->C")
                sim.new_reaction_counter(1,P_and_E_to_C,"P+E<-C")
            '''
            elif ii == "distribution":
                combistamp=self.stamp()
                distribution_sx= self.path + combistamp + "distributionSx.h5"
                distribution_sy= self.path + combistamp + "distributionSy.h5"
                distribution_sz=self.path + combistamp + "distributionSz.h5"
                distribution_px= self.path + combistamp + "distributionSx.h5"
                distribution_py= self.path + combistamp + "distributionSy.h5"
                distribution_pz=self.path + combistamp + "distributionSz.h5"
                number_S=self.path + combistamp + "number_s.h5"
                number_P=self.path + combistamp + "number_p.h5"
                number_C=self.path + combistamp + "number_c.h5"
                sim.new_particle_numbers(1,number_S,0)
                sim.new_particle_numbers(1,number_P,1)
                sim.new_particle_numbers(1,number_C,3)
                sim.new_probability_density(self.length, distribution_sx, 0, range(-self.boxsize/2,self.boxsize/2,50),0)
                sim.new_probability_density(self.length, distribution_sy, 0, range(-self.boxsize/2,self.boxsize/2,50),1)
                sim.new_probability_density(self.length, distribution_sz, 0, range(-self.boxsize/2,self.boxsize/2,50),2)

                sim.new_(range(100,200), radial_S, np.arange(self.reactiondistance, self.boxsize/2, (self.boxsize/2.0-self.reactiondistance)/50.0 ), np.array([[2,0],[3,0]]))
                sim.new_radial_distribution(range(100,200), radial_P, np.arange(self.reactiondistance, self.boxsize/2, (self.boxsize/2.0-self.reactiondistance)/50.0 ), np.array([[2, 1],[3,1]]))
                #sim.new_radial_distribution(self.length, radial_S,r_radial, np.array([[2,0],[3,0]]))
                #sim.new_radial_distribution(self.length, radial_P,r_radial, np.array([[2, 1],[3,1]]))
            '''

        sim.run(self.length)
        sim.write_observables_to_file()
        sim.delete_all_observables()

        #sim.show_world()

"""
particlenumber=1000
boxsize=7.0
D=1.0/6.0
R=1.0
alpha=0.5
lambda_plus=1.0
lambda_c=1.0
lambda_minus=1.0
tau=boxsize**(3)
tau= tau/(lambda_plus*(R*np.pi*particlenumber*6*D)**alpha)
print 2**13

complex=Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,alpha,1024,boxsize,particlenumber,0.1)
complex.observables=["all"]
complex.info()
complex.run()
"""