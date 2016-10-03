import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import lambertw

'''
def func(self,params, X, Y):
                    # extract current values of fit parameters from input array
                    K_m = params[0]
                    # compute chi-square
                    chi2 = 0.0
                    S_0=(self.particles*1.0)/(1.0*boxsize**3)
                    E_0=1.0/(1.0*boxsize**3)
                    k2=self.micro_reactionrate_complex
                    for n in range(len(X)):
                        x = X[n]
                        # The function y(x)=a+b*x+c*x^2 is a polynomial
                        # in this example.
                        y=K_m*lambertw((S_0/K_m)*np.exp((-k2*E_0*x+S_0)/K_m))
                        #S_t=S_t.astype(float)
                        chi2 = chi2 + (Y[n] - y)**2

'''
r=np.linspace(0,3.5)
v=7**3/50
dr=-r+(r**3+((v*3)/(4*np.pi)))**(1.0/3.0)
plt.plot(dr,r)
plt.show()

r_momentan=1
r_radial=[r_momentan]
for _ in range(30):
    dr=-r_momentan+(r_momentan**3+((v*3)/(4*np.pi)))**(1.0/3.0)
    r_momentan=r_momentan+dr
    r_radial.append(r_momentan)

plt.plot(r_radial,linestyle=" ",marker="+")
plt.show()

