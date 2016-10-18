import comlex_anal
import complex_sim
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw


particlenumber=200
boxsize=7.0
D=1.0/6.0
R=1.0
alpha=0.5
lambda_plus=1.0
lambda_c=1.0
lambda_minus=1.0
length=2**10


for ii in [0.06,0.11,0.81,1.41]:
    lambda_plus = ii
    lambda_minus= ii
    lambda_c =ii
    boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/200)+(4.0/3.0*np.pi))**(1/3.0)
    #complex=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.5,length,boxsize,particlenumber,0.05)
    #complex.observables=["number_of_particle","radial","MSD","reaction"]
    #for index1 in range(20):
    #    print R
    #    print index1
    #    complex.run()


    for index2 in range(400):
        complex1=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,1.0,length,boxsize,particlenumber,0.05)
        complex1.observables=["all"]
        print index2
        complex1.run()


for ii in [0.06,0.11,0.81,1.41]:
    lambda_plus = ii
    lambda_minus= ii
    lambda_c =ii
    boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/200)+(4.0/3.0*np.pi))**(1/3.0)
    #complexa=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.5,length,boxsize,particlenumber,0.05)
    #complexa.observables=["number_of_particle"]
    #complexa.show_observable()

    complexb=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,1.0,length,boxsize,particlenumber,0.05)
    complexb.observables=["number_of_particle"]
    complexb.show_observable()


    # Paramter
    boxsize=complexb.boxsize*1.0
    S_0=(complexb.particles*1.0)/(boxsize**3)
    E_0=1.0/(boxsize**3)
    k_minus=complexb.k_minus
    k_1=complexb.k_plus
    k_complex=complexb.k_complex
    K_m=(k_minus+k_complex)/k_1


    print k_1, K_m,k_minus,k_complex







    #t= np.linspace(0,length*complexb.timestep,length)


    S_t=K_m*lambertw((S_0/K_m)*np.exp((-k_complex*E_0*complexb.time+S_0)/K_m))
    ES_t=((E_0*S_t)/(K_m+S_t))*(1-np.exp(-(K_m+S_t)*k_1*complexb.time))
    E_t=E_0-ES_t
    P_t=S_0-S_t-ES_t

    ######

    plt.plot(complexb.time,ES_t/E_0,linestyle=":") #complex
    #plt.plot(t,E_t/E_0)  #enzyme
    plt.plot(complexb.time,S_t/S_0,linestyle=":")  #substrate
    plt.plot(complexb.time,P_t/S_0,linestyle=":")  #product


######


plt.ylabel('scaled concentration',fontsize=10)
plt.xlabel('scaled time',fontsize=10)




######






plt.legend(loc=4,fontsize='small')





plt.show()