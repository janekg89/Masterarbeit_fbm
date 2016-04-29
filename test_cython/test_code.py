import numpy
import genereatefracincrements as ginc
import matplotlib.pyplot as plt

__author__ = 'janek'
i = 1

N=2024

a = ginc.generateIncrements(N=N, D=numpy.array(0.3), tau=numpy.array(1.0), alpha=numpy.array(0.9), particles=3)
#for ii in range(3):
print a[0:100]

#plt.plot(a)
#plt.show()
'''
r_t=numpy.zeros(N)
r_t[1:]=numpy.cumsum(a)[:-1]
#plt.plot(a)
a = ginc.generateIncrements(N=N, D=numpy.array(1.3), tau=numpy.array(0.4), alpha=numpy.array(0.5))
a=numpy.array(a)
r_t=numpy.zeros(N)
print a.std()
r_t[1:]=numpy.cumsum(r_t)[:-1]
plt.plot(a,"g")
v_t=(numpy.random.normal(0,numpy.sqrt(0.4),size=N))
print v_t.std()
print numpy.sqrt(0.4)
plt.plot(v_t,"y")
plt.show()
'''






