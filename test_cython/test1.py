import numpy as np
import genereatefracincrements as ci
import matplotlib.pyplot as plt # module for plotting "a la" matlab
from mpl_toolkits.mplot3d import Axes3D

inc = ci.pyIncrements(500,1000)
inc.generateIncrements(0.5,0.5,0.5)
a = inc.returnIncrements()
print a.shape
a1d=a[:,:,0]
print a1d.shape



fig = plt.figure()



r_t=np.cumsum(a1d,axis=0)
r_t[0,:]=0
'''
#plt.plot(r_t[0,:])
#plt.plot(a[:,0,0])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(r_t[:,0], r_t[:,1], r_t[:,2])



print a1d.shape
'''


a1d1=a[:,:,1]


r_t1=np.cumsum(a1d1,axis=0)
r_t1[0,:]=0
#plt.plot(r_t[0,:])
#plt.plot(a[:,0,0])
ay = fig.add_subplot(111, projection='3d')
ay.plot(r_t1[:,0], r_t1[:,1], r_t1[:,2])


plt.show()
'''


plt.hist(a[:,0,0], 50, normed=1, facecolor='green', alpha=0.75)
plt.hist(a[0,:,0], 50, normed=1, facecolor='red', alpha=0.75)
plt.hist(a[0,0,:], 50, normed=1, facecolor='blue', alpha=0.75)




# add a 'best fit' line


plt.xlabel('Smarts')
plt.ylabel('Probability')
plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
plt.axis([40, 160, 0, 0.03])
plt.grid(True)

plt.show()

'''