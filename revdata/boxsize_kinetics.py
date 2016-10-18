import comlex_anal
import complex_sim
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw




D=1.0/6.0
R=1.0
alpha=0.5
lambda_plus=0.8
lambda_c=0.8
lambda_minus=0.8
length=2**10


for particlenumber in [50]:
    boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/200)+(4.0/3.0*np.pi))**(1/3.0)
    print boxsize

    #for index1 in range(500):
        #complex=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.5,length,boxsize,particlenumber,0.05)
        #complex.observables=["all"]
        #print R
        #print index1
        #complex.run()


    #for index2 in range(500):
     #   complex1=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,1.0,length,boxsize,particlenumber,0.05)
     #   complex1.observables=["number_of_particle","radial","reaction","MSD"]
     #   print index2
     #   complex1.run()






for particlenumber in [50]:
    #length=2**15
    boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/200)+(4.0/3.0*np.pi))**(1/3.0)
    complexa=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.5,length,boxsize,particlenumber,0.05)
    #complexa.observables=["number_of_particle"]
    #complexa.show_observable()
    #complexa.observables=["reaction"]
    #complexa.show_observable()

    complexa.observables=["radial"]
    complexa.show_observable()

    '''



    # Paramter
    boxsize=complexa.boxsize*1.0
    S_0=(complexa.particles*1.0)/(boxsize**3)
    E_0=1.0/(boxsize**3)
    k2=complexa.k_minus
    k_1=complexa.k_plus
    k_complex=complexa.k_complex
    K_m=(k2+k_complex)/k_1


    print k_1, K_m,k2,k_complex







    #t= np.linspace(0,length*complexb.timestep,length)


    S_t=K_m*lambertw((S_0/K_m)*np.exp((-k2*E_0*complexa.time+S_0)/K_m))
    ES_t=((E_0*S_t)/(K_m+S_t))*(1-np.exp(-(K_m+S_t)*k_1*complexa.time))
    E_t=E_0-ES_t
    P_t=S_0-S_t-ES_t

    ######

    plt.plot(complexa.time,ES_t/E_0,linestyle=":" ) #complex
    #plt.plot(t,E_t/E_0,linestyle=":")  #enzyme
    plt.plot(complexa.time,S_t/S_0,linestyle=":")  #substrate
    plt.plot(complexa.time,P_t/S_0,linestyle=":")  #product


    complexb=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.5,length,boxsize,particlenumber,0.05)
    complexb.observables=["number_of_particle"]
    complexb.show_observable()





    # Paramter
    boxsize=complexb.boxsize*1.0
    S_0=(complexb.particles*1.0)/(boxsize**3)
    E_0=1.0/(boxsize**3)
    k2=complexb.k_minus
    k_1=complexb.k_plus
    k_complex=complexb.k_complex
    K_m=(k2+k_complex)/k_1


    print k_1, K_m,k2,k_complex







    #t= np.linspace(0,length*complexb.timestep,length)


    S_t=K_m*lambertw((S_0/K_m)*np.exp((-k2*E_0*complexb.time+S_0)/K_m))
    ES_t=((E_0*S_t)/(K_m+S_t))*(1-np.exp(-(K_m+S_t)*k_1*complexb.time))
    E_t=E_0-ES_t
    P_t=S_0-S_t-ES_t

    ######

    plt.plot(complexb.time,ES_t/E_0,linestyle=":" ) #complex
    #plt.plot(t,E_t/E_0,linestyle=":")  #enzyme
    plt.plot(complexb.time,S_t/S_0,linestyle=":")  #substrate
    plt.plot(complexb.time,P_t/S_0,linestyle=":")  #product
    '''

    ######

plt.ylabel('radial distribution function',fontsize=10)
plt.xlabel('r',fontsize=10)
plt.legend(loc=4,fontsize='small')
#plt.savefig('../finalreport/data/checking_stability',dpi=300)
plt.show()

