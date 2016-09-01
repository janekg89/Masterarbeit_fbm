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


class RevreaDDy_analyse():
    def __init__(self, boxsize, Diffusion_E, Diffusion_S, Diffusion_P, reactiondistance, intrinsikreaction, length,
                 particles,alpha):
        self.boxsize = boxsize
        self.Diffusion_E = Diffusion_E
        self.Diffusion_S = Diffusion_S
        self.Diffusion_P = Diffusion_P
        self.reactiondistance = reactiondistance
        self.intrinsikreaction = intrinsikreaction
        self.length = length
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

    def simulation2(self):
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

    def simulation3(self):
        """

        :return: radialdistribution brownian
        """
        sim = rdy.Sim()
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

    def simulation4(self):
        """

        :return: number of substrate for brownian motion
        """
        sim = rdy.Sim()
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
        for _ in range(1):
            sim.add_particle(self.posenzyme, 2)
        sim.show_world()
        sim.delete_all_observables()
        numbers_a_name = "revdata/" + self.stamp() + "numbers_a.h5"
        sim.new_particle_numbers(1, numbers_a_name, 0)
        numbers_b_name = "revdata/" + self.stamp() + "numbers_b.h5"
        sim.new_particle_numbers(1, numbers_b_name, 1)
        print "NUMBERS A:", numbers_a_name
        print "NUMBERS B:", numbers_b_name
        sim.run(self.length)
        sim.write_observables_to_file()
        sim.delete_all_observables()
        return numbers_a_name, numbers_b_name

    def simulation5(self):
        """

        :return: number of substrates for fbm
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
        for _ in range(1):
            sim.add_particle(self.posenzyme, 2)
        sim.show_world()
        sim.delete_all_observables()
        numbers_a_name = "revdata/" + self.stamp() + "numbers_a.h5"
        sim.new_particle_numbers(1, numbers_a_name, 0)
        numbers_b_name = "revdata/" + self.stamp() + "numbers_b.h5"
        sim.new_particle_numbers(1, numbers_b_name, 1)
        print "NUMBERS A:", numbers_a_name
        print "NUMBERS B:", numbers_b_name
        sim.run(self.length)
        sim.write_observables_to_file()
        sim.delete_all_observables()
        fa = h5.File(numbers_a_name, 'r')
        fb = h5.File(numbers_b_name, 'r')
        return numbers_a_name, numbers_b_name
    def simulationtrajectory(self):
        sim = rdy.Sim()
        sim.kt = 2.437 # kJ/mol -> ~300 K
        sim.alpha = self.alpha
        sim.timestep = 1.0# ns
        sim.boxsize = self.boxsize# nm
        sim.delete_all_particle_types() # also deletes all interactions etc.
        sim.new_type("S", 0.01, self.Diffusion_S)
        sim.new_type("P", 0.01, self.Diffusion_P)
        sim.new_type("E", 0.01, self.Diffusion_E)
        sim.delete_all_reactions()
        sim.new_enzymatic("E+S<->E+S", 0, 1, 2, self.intrinsikreaction, 0.1, self.reactiondistance)
        sim.delete_all_particles()
        for i in range(self.particles):
            pos = self.positions[i,:]
            sim.add_particle(pos, 0)
            #pos = np.random.random(3) * 9.9 - 4.95
            #sim.add_particle(pos, 1)
            #pos = np.random.random(3) * 9.9 - 4.95
            #sim.add_particle(pos, 2)
        for _ in range(1):
            sim.add_particle(self.posenzyme, 2)
        sim.show_world()
        sim.delete_all_observables()
        mean_e_name = "revdata/" +self.stamp()+"traject.h5"
        sim.new_trajectory_unique(1,self.length,mean_e_name)
        sim.run(self.length)
        sim.write_observables_to_file()
        sim.delete_all_observables()
        fc = h5.File(mean_e_name, 'r')
        return  mean_e_name
    def simulationtrajectory2(self):
        sim = rdy.Sim("FractionalDiffusion")
        sim.kt = 2.437 # kJ/mol -> ~300 K
        sim.alpha = self.alpha
        sim.timestep = 1.0# ns
        sim.boxsize = self.boxsize# nm
        sim.delete_all_particle_types() # also deletes all interactions etc.
        sim.new_type("S", 0.01, self.Diffusion_S)
        sim.new_type("P", 0.01, self.Diffusion_P)
        sim.new_type("E", 0.01, self.Diffusion_E)
        sim.delete_all_reactions()
        sim.new_enzymatic("E+S<->E+S", 0, 1, 2, self.intrinsikreaction, 0.1, self.reactiondistance)
        sim.delete_all_particles()
        for i in range(self.particles):
            pos = self.positions[i,:]
            sim.add_particle(pos, 0)
            #pos = np.random.random(3) * 9.9 - 4.95
            #sim.add_particle(pos, 1)
            #pos = np.random.random(3) * 9.9 - 4.95
            #sim.add_particle(pos, 2)
        for _ in range(1):
            sim.add_particle(self.posenzyme, 2)
        sim.show_world()
        sim.delete_all_observables()
        mean_e_name = "revdata/" +self.stamp()+"traject.h5"
        sim.new_trajectory_unique(1,self.length,mean_e_name)
        sim.run(self.length)
        sim.write_observables_to_file()
        sim.delete_all_observables()
        return  mean_e_name


    def simulationreactionamount(self):
        """

        :return: number of substrate for brownian motion
        """
        sim = rdy.Sim()
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
        for _ in range(1):
            sim.add_particle(self.posenzyme, 2)
        sim.show_world()
        sim.delete_all_observables()
        numbers_a_name = "revdata/" + self.stamp() + "numberreaction_a.h5"
        sim.new_reaction_counter(1, numbers_a_name,"E+S<->E+S")
        sim.run(self.length)
        sim.write_observables_to_file()
        sim.delete_all_observables()
        return numbers_a_name


    def simulationreactionamountfrac(self):
        """

        :return: number of substrates for fbm
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
        for _ in range(1):
            sim.add_particle(self.posenzyme, 2)
        sim.show_world()
        sim.delete_all_observables()
        numbers_a_name = "revdata/" + self.stamp() + "numberreaction_a.h5"
        sim.new_reaction_counter(1, numbers_a_name,"E+S<->E+S")
        sim.run(self.length)
        sim.write_observables_to_file()
        sim.delete_all_observables()
        particles=sim.show_world()
        return numbers_a_name, particles
    def simulationpositions(self):
        """

        :return: number of substrates for fbm
        """
        sim = rdy.Sim("FractionalDiffusion")
        sim.kt = 2.437  # kJ/mol -> ~300 K
        sim.alpha = 0.5
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
        for _ in range(1):
            sim.add_particle(self.posenzyme, 2)
        sim.show_world()
        sim.delete_all_observables()
        numbers_a_name = "revdata/" + self.stamp() + "numberreaction_a.h5"
        sim.new_reaction_counter(1, numbers_a_name)
        sim.run(self.length)
        num_particles = sim.world.getNumberOfParticles()
        result = []
        result2 = []
        for i in range(num_particles):
            if sim.world.getTypeId(i) == 0:
                pos = sim.world.getPosition(i)
                result.append(pos)
            elif sim.world.getTypeId(i) == 1:
                pos = sim.world.getPosition(i)
                result2.append(pos)
        return  np.array(result),np.array(result2)

Rev = RevreaDDy_analyse(30, 0.2, 0.2, 0.2, 2.0, 0.5, 16384/4, 500,0.5)




"""
fd=Rev.simulation2()
fc=Rev.simulation3()


#,fc
plt.plot(fd['binCenters'],fd['bins'], label="fractional")
plt.plot(fc['binCenters'],fc['bins'], label="brownian")

plt.ylabel('Bins', fontsize=10)
plt.xlabel('Bin-Centers', fontsize=10)
plt.legend()
#plt.savefig('finalreport/data/fractionalradialdistributionrevreaddy.png',dpi=300)   # save the figure to file
plt.show()


-----------------

colors = iter(cm.rainbow(np.linspace(0, 1, 4)))
colornow = next(colors)
falist = []
fclist = []
faflist = []
names_a=[]
names_b=[]
names_af=[]
names_bf=[]
for ii in range(10):
    print ii
    numbers_a_name, numbers_b_name = Rev.simulation4()  # ,fc
    numbers_a_namef, numbers_b_namef= Rev.simulation5()  # ,fc
    names_a.append(numbers_a_name)
    names_b.append(numbers_b_name)
    names_af.append(numbers_a_namef)
    names_bf.append(numbers_b_namef)
for i in range(10) :
    fa = h5.File(names_a[i], 'r')
    fb = h5.File(names_b[i], 'r')
    faf = h5.File(names_af[i], 'r')
    fbf = h5.File(names_bf[i], 'r')
    falist.append(fa['particleNumbers'])
    faflist.append(faf['particleNumbers'])
falist = np.array(falist)
faflist = np.array(faflist)
paertnum = falist.mean(axis=0)
paertnumf = faflist.mean(axis=0)
plt.plot(fbf['times'], paertnumf, color=colornow, label="$D_S+D_E$=$%.2f$ fractional" % (0.4))
plt.plot(fb['times'], paertnum, color=next(colors), label="$D_S+D_E$=$%.2f$ brownian" % (0.4))
plt.ylabel('S(t)', fontsize=10)
plt.xlabel('t in ns', fontsize=10)
plt.legend()
plt.show()
------------------------

traj=Rev.simulationtrajectory()

#traj['positions']
#fb=Rev.simulationtrajectory2()
with open('revdata/myfile.xyz','w') as f:
    for length in range(traj['positions'].shape[0]):
        s=str(traj['positions'].shape[1])+"\n"
        f.write(s)
        f.write('comment\n')
        for particle in range(traj['positions'].shape[1]):
            st="S"+"          "+str(traj['positions'][length,particle,0])+"          "+str(traj['positions'][length,particle,1])+"          "+str(traj['positions'][length,particle,2])+'\n'
            f.write(st)
            # python will convert \n to os.linesep
----------------------------------------

change alpha
"""
colors = iter(cm.rainbow(np.linspace(0, 1, 20)))

for iii in [0.5,0.6,0.7,0.8,0.9]:
    Rev.alpha=iii
    colornow = next(colors)
    fclist = []
    faflist = []
    names_a=[]
    names_b=[]
    names_af=[]
    names_bf=[]
    for ii in range(2):
        print ii
        numbers_a_namef, particles= Rev.simulationreactionamountfrac()  # ,fc
        names_af.append(numbers_a_namef)
    for i in range(2):
        faf = h5.File(names_af[i], 'r')
        faflist.append(faf['forwardCounter'][:])
    faflist = np.array(faflist)
    reacnumf = faflist.mean(axis=0)
    plt.plot(faf['times'],reacnumf, color= next(colors), label="$alpha$=$%.2f$ fractional" % (iii))
    plt.ylabel('Number of forward reactions', fontsize=10)
    plt.xlabel('t in ns', fontsize=10)
    plt.legend()
names_a=[]
falist = []
for ii in range(2):
        print ii
        numbers_a_name = Rev.simulationreactionamount()  # ,fc
        names_a.append(numbers_a_name)
for i in range(2):
        fa = h5.File(names_a[i], 'r')
        falist.append(fa['forwardCounter'][:])
falist = np.array(falist)
reacnum = falist.mean(axis=0)
plt.plot(fa['times'], reacnum, color=next(colors),label="$alpha$=$%.2f$ brownian" % (iii))
plt.show()
"""
----------------------------------

num_particles,num_particles2=Rev.simulationpositions()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(num_particles[:,0], num_particles[:,1],num_particles[:,2], c='b')
ax.scatter(num_particles2[:,0], num_particles2[:,1], num_particles2[:,2], c='r')
plt.show()

"""
