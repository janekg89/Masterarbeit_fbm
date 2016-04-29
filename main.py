# -*- coding: utf-8 -*-
__author__ = 'janekg89'
import numpy as np # module for scientific computing
import matplotlib.pyplot as plt # module for plotting "a la" matlab
import simulation
import timeit
import analyse_tool
import test_cython.genereatefracincrements as ginc

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
    length=1000
    steps=30
    gaussianparamter=[]
    for i in range(steps):
        c=analyse_tool.Analyse(D=2,particles=5000,length=length+1,alpha=0.5,dt=1)
        shades=np.linspace(0,1,steps+1)
        b=c.nongaussian_parameter()
        gaussianparamter.append(b)
    gaussianparamter=np.array(gaussianparamter)
    plt.errorbar(range(length+1), gaussianparamter.mean(axis=0), yerr=gaussianparamter.std(axis=0))
    plt.xlabel('t', fontsize=14)
    plt.ylabel('$\\alpha_2(t)$ non-Gaussian parameter', fontsize=14)
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
        #d=c.rescaled_analytical_distribution(t=j,r_dis=50)
        h=c.rescaled_function(t=j,histpoints=35)
        plt.loglog(h[0],h[1],color='%f' %(shade), label="rescaled function at $t=%d$" %(j))
        #plt.plot(d[0],d[1],color='%f' %(shade))
    plt.xlabel('rescaled distance $ r_{res} $ ', fontsize=14)
    plt.ylabel('rescaled distribution', fontsize=14)
    plt.legend(loc=1)
    plt.show()

def show_distrib():

    colornum=-1
    length=1000
    steps=50
    shades=np.linspace(0,1,steps+1)
    c=analyse_tool.Analyse(D=2,particles=6000,length=length+1,alpha=0.5,dt=1,)
    for j in range(100,1000,900/steps):
        colornum=1+colornum
        shade=shades[colornum]
        f=c.analytical_distribution_of_particles(t=j,r_dis=50)
        g=c.distribution(t=j)
        plt.semilogy(g[0],g[1],color='%f' %(shade))
        plt.plot(f[0],f[2],color='%f' %(shade))
    plt.show()

def plot_msd_invert():
    c=analyse_tool.Analyse(D=2,particles=500,length=5000,alpha=0.5,dt=1)
    c.invert_time()
    plt.show()
#c=analyse_tool.Analyse(D=2,particles=6000,length=1000,alpha=0.5,dt=1)
#msddopp=[]
#for i in range(5999):
#    msdd=c.msd_time(i=i)
#    msddopp.append(msdd)
#msddopp=np.array(msddopp)
#c.trajectory=msddopp
#c.plotting(showlegend="yes")
#c.plotting(showlegend="yes")
#plt.show()
#plot_msd_invert()
#show_gaussian()
#show_rescaled()
#show_distrib()
#show_rescaled()
#def plot_different_dt():
#        gaussianparamter=[]
#        for i in range(steps):
 #           c=analyse_tool.Analyse(D=2,particles=5000,length=length+1,alpha=0.5,dt=j)
  #          shades=np.linspace(0,1,steps+1)
   #         b=c.nongaussian_parameter()
    #        gaussianparamter.append(b)
     #   gaussianparamter=np.array(gaussianparamter)
      #  plt.errorbar(range(length+1), gaussianparamter.mean(axis=0), yerr=gaussianparamter.std(axis=0))
    #plt.show()
#plot_different_dt()



def plot_ensemble_mean_of_time_msd():
        c=analyse_tool.Analyse(D=2,particles=10,length=1000,alpha=0.5,dt=1,version="cpp")
        msd_all=[]
        for ii in range(c.particles):
            print ii
            msd,std=c.msd_time(i=ii)
            msd_all.append(msd)
        colors=['r','b','g','k','c','w','b','r','g','b','k','c','w','b','r','g','b','k','c','w','bo','ro','go','bo','ko','co','wo','bo']
        msd_time=np.array(msd_all)
        msd_time_mean=msd_time.mean(axis=0)
        plt.plot(c.t*c.dt,c.msdanalyt(),":",color=colors[1], label="analytisch D=%f,particles=%d,length=%d,alpha=%f" %(c.D,c.particles,c.n,c.alpha))
        #plt.errorbar(c.t*c.dt, msd_time_mean, yerr=0,label="Spektrale Methode mit D=%f,particles=%d, length=%d ,alpha=%f" %(c.D,c.particles,c.n,c.alpha))
        #plt.show()
        return msd_time_mean



#plot_msd_invert()


#plt.plot(plot_ensemble_mean_of_time_msd())
#plt.show()
#show_gaussian()

#show_rescaled()


#show_distrib()


print timeit.timeit("import analyse_tool; analyse_tool.Analyse(D=2,particles=1024,length=1024,alpha=0.5,dt=1,version='cpp').compute_trajectory()", number=1)
print timeit.timeit("import analyse_tool; analyse_tool.Analyse(D=2,particles=1024,length=1024,alpha=0.5,dt=1).compute_trajectory()", number=1)

#e=analyse_tool.Analyse(D=2,particles=1000,length=1024*2,alpha=0.5,dt=1.,version='cpp')
#plt.plot(e.compute_trajectory()[0],"r")
'''
a = ginc.generateIncrements(N=e.n, D=np.array(2.0), tau=np.array(1.0), alpha=np.array(0.5))

v_t=(np.random.normal(0,np.sqrt(1.),size=e.n))
v_t=np.array(v_t)
v_frq=np.fft.fft(v_t)
v_ano_frq= np.sqrt(e.z().real[:2048]*2.)*v_frq
plt.plot(v_ano_frq.real)
plt.plot(a,"r")
print v_ano_frq.std()
print np.array(a).std()


d=analyse_tool.Analyse(D=2.0,particles=2000,length=1000,alpha=0.5,dt=0.5,x=4)
 #plt.plot(d.compute_trajectory()[0])
e=analyse_tool.Analyse(D=2,particles=2000,length=1000,alpha=0.5,dt=1.,version='cpp')
print d.compute_trajectory().shape
print e.compute_trajectory().shape

d.plotting(showlegend="Yes")
#e.plotting(msdtype="time",scale="lin",showlegend="Yes")
e.plotting(showlegend="Yes",scale="lin")
#d.plotting(scale="lin",showlegend="Yes")
plt.show()

import revreaddy.build.Debug.revreaddy as rdy
s =rdy.Sim()
s.run(steps=10000, timestep=0.1)
'''