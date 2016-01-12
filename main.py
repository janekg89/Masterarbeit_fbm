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
length=1000
steps=20
c=analyse_tool.Analyse(D=2,particles=9000,length=length+1,alpha=0.5)
shades=np.linspace(0,1,steps+1)
colornum=-1
for j in range(900,1000,50): #TODO find out why to lines if you run it
    print j
    colornum=1+colornum
#    d=c.rescaled_analytical_distribution(t=j,r_dis=50)
    e=c.rescaled_analytical_distribution(t=j,r_dis=50)
    #f=c.analytical_distribution_of_particles(t=j,r_dis=50)
    #g=c.distribution(t=j)
    h=c.rescaled_function(t=j,histpoints=35)
    shade=shades[colornum]
    print shade
    #f=c.rescaled_function(j)
    #plt.plot(g[0],g[1],'b')
    #plt.plot(f[0],f[2],'b')
    plt.plot(e[0],e[1])
    #plt.plot(f[0],f[1],'g')
    plt.plot(h[0],h[1],color='%f' %(shade))
#plt.plot(b,a)
#plt.plot(c[0],c[1])
#plt.show()
#a,b=analyse_tool.Analyse(D=2,particles=200,length=2000,alpha=0.4).distribution(t=1)
#plt.plot(b,a)
plt.show()
#print b.size




