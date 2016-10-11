import comlex_anal
import complex_sim
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw





particlenumber=20
boxsize=8.0
D=1.0/6.0
R=1.0
alpha=1.0
lambda_plus=1.0
lambda_c=1.0
lambda_minus=1.0
length=2**14
tau=0.05

def show_particle(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau):
        complexb=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
        complexb.observables=["number_of_particle"]
        complexb.show_observable()

        boxsize1=complexb.boxsize*1.0
        S_0=(complexb.particles*1.0)/(boxsize1**3)
        E_0=1.0/(boxsize1**3)
        k2=complexb.k_minus
        k_1=complexb.k_plus#*(complexb.time*0.05)**(np.log(complexb.alpha)/np.pi)
        k_complex=complexb.k_complex
        K_m=(k2+k_complex)/k_1
        k_erban=4.0*np.pi*complexb.Diffusion*(complexb.reactiondistance-(np.sqrt(complexb.Diffusion/complexb.micro_reactionrate_forward)*np.tanh(complexb.reactiondistance*np.sqrt(complexb.micro_reactionrate_forward/complexb.Diffusion))))
        print "erban chapmann:k1=",k_erban


        S_t=K_m*lambertw((S_0/K_m)*np.exp((-k_complex*E_0*complexb.time+S_0)/K_m))
        ES_t=((E_0*S_t)/(K_m+S_t))*(1-np.exp(-(K_m+S_t)*k_1*complexb.time))
        E_t=E_0-ES_t
        P_t=S_0-S_t-ES_t


        #plt.plot(complexb.time,ES_t/E_0,linestyle=":" ) #complex
        #plt.plot(t,E_t/E_0,linestyle=":")  #enzyme
        plt.plot(complexb.time,S_t/S_0,linestyle=":")  #substrate
        #plt.plot(complexb.time,P_t/S_0,linestyle=":")  #product

        #print k_1, K_m,k2,k_complex

        #plt.plot(complexb.var_all)

        #plt.errorbar(range(len(complexb.k_all)),complexb.k_all,complexb.var_all)

def show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau):
        complexb=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
        complexb.observables=["radial"]
        complexb.show_observable()
def show_reaction(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau):
        complexb=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
        complexb.observables=["reaction"]
        complexb.show_observable()
def show_MSD(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau):
        complexb=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
        complexb.observables=["MSD"]
        complexb.show_observable()

        plt.plot(complexb.time[1:],6*((complexb.time[1:])**complexb.alpha)*complexb.Diffusion,linestyle="--", label="analytical")
"""
for boxsize in [8,9,10,11]:
    alpha=1.0
    show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
    #alpha=1.0
    #show_particle(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)

for lambda_plus in [1.0]:
    alpha=1.0
    #show_particle(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
    show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
    alpha=0.5
    #show_particle(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
    show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)

    #alpha=0.5
    #show_particle(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)

for lambda_plus in [1.00]:
    alpha=1.0
    show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
    #alpha=0.5
    #show_particle(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
"""
for lambda_minus in [1.00]:
     alpha=0.5
     show_reaction(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
"""
for lambda_c in [0.1,1.00]:
     alpha=0.5
     show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
     #alpha=0.5
     #show_particle(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)



for boxsize in [10.5]:
    #alpha=0.5
    #show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
    alpha=1.0
    show_particle(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
"""
plt.ylabel('radial distribution function',fontsize=10)
plt.xlabel('r',fontsize=10)
plt.legend(loc=4,fontsize='small')
#plt.savefig('../finalreport/data/checking_stability',dpi=300)
plt.show()


