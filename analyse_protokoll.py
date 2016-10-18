# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import comlex_anal
import complex_sim
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw

"""
Standart Starting Parameter
"""
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
        """

        :param D:
        :param R:
        :param lambda_plus:
        :param lambda_minus:
        :param lambda_c:
        :param alpha:
        :param length:
        :param boxsize:
        :param particlenumber:
        :param tau:
        :return: plots scaled concentration over time with a fitt of PDE od Michaelis Menten
        """
        complexb=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
        complexb.observables=["number_of_particle"]
        complexb.show_observable()

        boxsize1=complexb.boxsize*1.0
        S_0=(complexb.particles*1.0)/(boxsize1**3)
        E_0=1.0/(boxsize1**3)
        k2=complexb.k_minus
        k_1=complexb.k_plus
        k_complex=complexb.k_complex
        K_m=(k2+k_complex)/k_1

        #k_erban=4.0*np.pi*complexb.Diffusion*(complexb.reactiondistance-(np.sqrt(complexb.Diffusion/complexb.micro_reactionrate_forward)*np.tanh(complexb.reactiondistance*np.sqrt(complexb.micro_reactionrate_forward/complexb.Diffusion))))
        #print "erban chapmann:k1=",k_erban


        S_t=K_m*lambertw((S_0/K_m)*np.exp((-k_complex*E_0*complexb.time+S_0)/K_m))
        ES_t=((E_0*S_t)/(K_m+S_t))*(1-np.exp(-(K_m+S_t)*k_1*complexb.time))
        E_t=E_0-ES_t
        P_t=S_0-S_t-ES_t


        plt.plot(complexb.time,ES_t/E_0,linestyle=":" ) #complex
        plt.plot(complexb.time,E_t/E_0,linestyle=":")  #enzyme
        plt.plot(complexb.time,S_t/S_0,linestyle=":")  #substrate
        plt.plot(complexb.time,P_t/S_0,linestyle=":")  #product



def show_k1(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau):
        """

        :param D:
        :param R:
        :param lambda_plus:
        :param lambda_minus:
        :param lambda_c:
        :param alpha:
        :param length:
        :param boxsize:
        :param particlenumber:
        :param tau:
        :return:
        """
        complexb=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
        
        complexb.observables=["number_of_particle"]
        complexb.show_observable(show_plot="off")

        boxsize1=complexb.boxsize*1.0
        S_0=(complexb.particles*1.0)/(boxsize1**3)
        E_0=1.0/(boxsize1**3)
        k2=complexb.k_minus
        k_1=complexb.k_plus#*(complexb.time*0.05)**(np.log(complexb.alpha)/np.pi)
        k_complex=complexb.k_complex
        K_m=(k2+k_complex)/k_1
        k_erban=4.0*np.pi*complexb.Diffusion*(complexb.reactiondistance-(np.sqrt(complexb.Diffusion/(complexb.micro_reactionrate_forward))*np.tanh(complexb.reactiondistance*np.sqrt((complexb.micro_reactionrate_forward)/complexb.Diffusion))))
        print "erban chapmann:k1=",k_erban
        print "fitted k1",k_1

def show_k2(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau):
        complexb=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)

        complexb.observables=["kminus"]
        complexb.show_observable()
def show_particle_alpha(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau,k1_alpha_0):
        complexb=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)

        complexb.observables=["Michaelis_fit_alpha"]

        complexb.show_observable(k1_alpha_0=k1_alpha_0)
        return complexb.k1_alpha_0
def show_k1_better(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau):
        complexb=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)

        complexb.observables=["kplus"]
        complexb.show_observable()

        k_erban=4.0*np.pi*complexb.Diffusion*(complexb.reactiondistance-(np.sqrt(complexb.Diffusion/(complexb.micro_reactionrate_forward))*np.tanh(complexb.reactiondistance*np.sqrt((complexb.micro_reactionrate_forward)/complexb.Diffusion))))
        print "erban chapmann:k1=",k_erban
        k_1=complexb.k_plus
        print "fitted k1",k_1


def show_kcomplex(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau):
        complexb=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)

        complexb.observables=["kcomplex"]
        complexb.show_observable()
       
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
for alpha in [0.5,0.6,0.7,0.8,0.9,1.0]:

    #show_reaction(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
    show_k1_better(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
    #show_particle_alpha(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
    

    #,,,alpha=1.0
plt.ylabel('k1',fontsize=10)
plt.xlabel('t',fontsize=10)
plt.legend(loc=4,fontsize='small')
#plt.savefig('../finalreport/data/checking_stability',dpi=300)
plt.show()

k1_alpha_0=0
for alpha in [1.0,0.9,0.5]:

    #show_reaction(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
    #show_k1_better(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
    print k1_alpha_0, "wichitger test"
    if alpha==1.0:
        k1=show_particle_alpha(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,8,particlenumber,tau,k1_alpha_0)
    else:
        k1=show_particle_alpha(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau,k1_alpha_0)

    if alpha == 1.0:
        k1_alpha_0=k1



plt.ylabel('scaled concentration',fontsize=10)
plt.xlabel('scaled time',fontsize=10)
plt.legend(loc=4,fontsize='small')
#plt.savefig('../finalreport/data/checking_stability',dpi=300)
plt.show()

for boxsize in [5,8.0,11]:

    #show_reaction(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
    show_k1_better(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)


    #,,,alpha=1.0
plt.ylabel('k1',fontsize=10)
plt.xlabel('t',fontsize=10)
plt.legend(loc=4,fontsize='small')
#plt.savefig('../finalreport/data/checking_stability',dpi=300)
plt.show()

   #show_particle(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)


for alpha in [0.5,0.6,0.7,0.8,0.9,1.0]:

    #show_reaction(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
    show_k2(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)



plt.ylabel('kminus',fontsize=10)
plt.xlabel('t',fontsize=10)
plt.legend(loc=4,fontsize='small')
#plt.savefig('../finalreport/data/checking_stability',dpi=300)
plt.show()





for alpha in [0.5,0.6,0.7,0.8,0.9,1.0]:

    #show_reaction(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
    show_kcomplex(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)



plt.ylabel('kcomplex',fontsize=10)
plt.xlabel('t',fontsize=10)
plt.legend(loc=4,fontsize='small')
#plt.savefig('../finalreport/data/checking_stability',dpi=300)
plt.show()




for boxsize in [8.0]:

    #show_reaction(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
    show_k2(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)



plt.ylabel('kminus',fontsize=10)
plt.xlabel('t',fontsize=10)
plt.legend(loc=4,fontsize='small')
#plt.savefig('../finalreport/data/checking_stability',dpi=300)
plt.show()


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

for lambda_minus in [1.00]:
     alpha=1.0
     show_reaction(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)

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



for boxsize in [8.001]:
   show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
plt.legend()
plt.show()

for boxsize in [14,15]:
    show_k1_better(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
plt.legend()
plt.show()
"""
for lambda_minus in [0.0]:
    lambda_plus=0.1
    show_particle_alpha(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau,0)
plt.legend()
plt.show()