import numpy
import genereatefracincrements as ginc

__author__ = 'janek'
N=2024
a = ginc.generateIncrements(N=N, D=numpy.array(0.3), tau=numpy.array(1.0), alpha=numpy.array(0.9), particles=3)
print a[0:100]







