import numpy
import genereatefracincrements as ginc

__author__ = 'janek'
a = ginc.generateIncrements(N=100, D=numpy.array(0.3), tau=numpy.array(0.3), alpha=numpy.array(0.5))

print a
