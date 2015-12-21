# -*- coding: utf-8 -*-
__author__ = 'janekg89'
import numpy as np # module for scientific computing
import matplotlib.pyplot as plt # module for plotting "a la" matlab
import simulation
import timeit
import analyse_tool

#print timeit.timeit(" import Spectral_method; Spectral_method.Felix_Method(D=2,particles=1,length=int(1e6),alpha=0.5).compute_trajectory()", number=1)
#Spectral_method.Felix_Method(D=2,particles=1,length=1000,alpha=0.5).compute_trajectory()
analyse_tool.Analyse(D=2,particles=1000,length=1000,alpha=0.6).plotting()
analyse_tool.Analyse(D=2,particles=1000,length=1000,alpha=0.7).plotting()
analyse_tool.Analyse(D=2,particles=1000,length=1000,alpha=1).plotting()
plt.show()
