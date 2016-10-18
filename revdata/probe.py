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
alpha=[1.0,0.9,0.8,0.7,0.6,0.5]
alpha=np.array(alpha)
h=[0.0,0.013,0.033,0.06,0.097,0.22]
k_0=[1.0,]
print (np.log(h[5])-0)/(-0.5)

plt.plot(h)
plt.plot(-0.05*np.exp(1)+0.05*np.exp(1/alpha))
plt.show()

