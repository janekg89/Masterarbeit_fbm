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
        plt.plot(range(2**14/20),np.ones(2**14/20)*k_erban,linestyle="--",label="Erban Chapmann")
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
for lambda_minus in [0.0]:
    lambda_plus=0.1
    #boxsize=15
    show_particle_alpha(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau,0)
plt.ylabel('scaled concentrations', fontsize=10)
plt.xlabel('scaled time', fontsize=10)
plt.legend(loc=2)
plt.savefig('./finalreport/data/chapman-limit-concentrations1',dpi=300)
plt.fill_between(np.linspace(100/20,18000/20),0.8*10**(-5), 10**(-5),color='blue', alpha=0.5)
plt.show()

for lambda_minus in [0.0]:
    lambda_plus=0.1
    #boxsize=15
    show_k1_better(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
plt.ylabel('$k_{+}$', fontsize=12)
plt.xlabel('scaled time', fontsize=10)
plt.legend(loc=3,fontsize='small')
plt.savefig('./finalreport/data/chapman-limit-concentrations1_k1',dpi=300)
#plt.k1_better(np.linspace(100/20,18000/20),0.8*10**(-5), 10**(-5),color='blue', alpha=0.5)
plt.show()


"""
for boxsize in [8.0]:
    lambda_minus=0.0
    lambda_plus=0.1
    show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
plt.ylabel('$r(t)$', fontsize=12)
plt.xlabel('scaled time', fontsize=10)
plt.legend(loc=3,fontsize='small')
plt.savefig('./finalreport/data/chapman-limit-radial',dpi=300)
plt.show()
