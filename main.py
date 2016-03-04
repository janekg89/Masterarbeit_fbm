# -*- coding: utf-8 -*-
__author__ = 'janekg89'
import numpy as np # module for scientific computing
import matplotlib.pyplot as plt # module for plotting "a la" matlab
import simulation
import timeit
import analyse_tool
import cProfile
import pstats
#print timeit.timeit(" import Spectral_method; Spectral_method.Felix_Method(D=2,particles=1,length=int(1e6),alpha=0.5).compute_trajectory()", number=1)
#Spectral_method.Felix_Method(D=2,particles=1,length=1000,alpha=0.5).compute_trajectory()
#analyse_tool.Analyse(D=2,particles=1000,length=1000,alpha=0.6).plotting()


#a.plotting(error=2)
#for j in range(10):
#    a.plotting(msdtype="time",particlemsdtime=j)
#a.plotting(error=2)
#plt.show()
#for j in range(10):
#a,b=analyse_tool.Analyse(D=2,particles=5000,length=200,alpha=0.5).distribution(t=199)
#print analyse_tool.Analyse(D=2,particles=50,length=20,alpha=0).D
def show_gaussian():
    length=100
    steps=30
    gaussianparamter=[]
    for i in range(steps):
        c=analyse_tool.Analyse(D=2,particles=5000,length=length+1,alpha=0.5,dt=1)
        shades=np.linspace(0,1,steps+1)
        b=c.nongaussian_parameter()
        gaussianparamter.append(b)
    gaussianparamter=np.array(gaussianparamter)
    plt.errorbar(range(length+1), gaussianparamter.mean(axis=0), yerr=gaussianparamter.std(axis=0))
    plt.show()

def show_rescaled():

    colornum=-1
    length=1100
    steps=5
    shades=np.linspace(0,1,steps+1)
    c=analyse_tool.Analyse(D=2,particles=10000,length=length+1,alpha=0.5,dt=1)
    for j in range(100,1100,1000/steps):
        colornum=1+colornum
        shade=shades[colornum]
        d=c.rescaled_analytical_distribution(t=j,r_dis=50)
        h=c.rescaled_function(t=j,histpoints=35)
        plt.plot(h[0],h[1],color='%f' %(shade), label='%f' %(shade))
        plt.plot(d[0],d[1],color='%f' %(shade))

#plt.legend(loc=1)
    plt.show()

def show_distrib():

    colornum=-1
    length=1000
    steps=50
    shades=np.linspace(0,1,steps+1)
    c=analyse_tool.Analyse(D=2,particles=6000,length=length+1,alpha=0.5,dt=1)
    c.invert_time()
    for j in range(100,1000,900/steps):
        colornum=1+colornum
        shade=shades[colornum]
        f=c.analytical_distribution_of_particles(t=j,r_dis=50)
        g=c.distribution(t=j)
        plt.semilogy(g[0],g[1],color='%f' %(shade))
        plt.plot(f[0],f[2],color='%f' %(shade))
    plt.show()

def plot_msd_invert():
    c=analyse_tool.Analyse(D=2,particles=6000,length=1000,alpha=0.5,dt=1)
    c.invert_time()
    c.plotting(showlegend="yes")
    c.plotting("time",showlegend="yes")
    plt.show()
c=analyse_tool.Analyse(D=2,particles=6000,length=1000,alpha=0.5,dt=1)
msddopp=[]
for i in range(5999):
    msdd=c.msd_time(i=i)
    msddopp.append(msdd)
msddopp=np.array(msddopp)
c.trajectory=msddopp
c.plotting(showlegend="yes")
c.plotting("time",showlegend="yes")
plt.show()
#plot_msd_invert()
#show_gaussian()
#show_rescaled()
#show_distrib()



