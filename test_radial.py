import comlex_anal
import complex_sim
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw
particlenumber=200
boxsize=3.5
D=1.0/6.0
R=1.0
length=2**5
lambda_plus=10000.0
lambda_c=1000.0
lambda_minus=0.0
boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/200)+(4.0/3.0*np.pi))**(1/3.0)
for index1 in range(100):
      print index1
      complex=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,1.0,length,boxsize,particlenumber,0.05)
      complex.observables=["radial"]
      complex.run()

complexa=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,1.0,length,boxsize,particlenumber,0.05)
complexa.observables=["number_of_particle"]
complexa.show_observable()
#lambda_c=0.1
#complexa=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,1.0,length,boxsize,particlenumber,0.05)
#complexa.observables=["radial"]
#complexa.show_observable()

plt.ylabel('$ \\frac{N}{N_{max}} $',fontsize=10)
plt.xlabel('t in $\\Gamma$',fontsize=10)
plt.legend(loc=4,fontsize='small')
#plt.savefig('../finalreport/data/checking_stability',dpi=300)
plt.show()
