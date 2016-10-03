import comlex_anal
import complex_sim
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw



boxsize=7.0
D=1.0/6.0
R=1.0
alpha=0.5
lambda_plus=1.0
lambda_c=1.0
lambda_minus=1.0
length=2**11


for particlenumber in [80]:
    lambda_plus=1
    lambda_c=0.1
    lambda_minus=0.1
    length=2**15
    boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/200)+(4.0/3.0*np.pi))**(1/3.0)
    print boxsize
    complex=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.5,length,boxsize,particlenumber,0.05)
    complex.observables=["number_of_particle","radial","MSD","reaction"]
    for index1 in range(500):
        print index1
        complex.run()



'''
for particlenumber in [10,50,100,500,1000]:
    boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/200)+(4.0/3.0*np.pi))**(1/3.0)
    print boxsize
    #complex=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.5,length,boxsize,particlenumber,0.05)
    #complex.observables=["number_of_particle","radial","MSD","reaction"]
    #for index1 in range(500):
    #    print R
    #    print index1
    #    complex.run()

    complex1=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,1.0,length,boxsize,particlenumber,0.05)
    complex1.observables=["number_of_particle","radial","reaction","MSD"]
    #for index2 in range(50):
     #   print index2
    #    complex1.run()

'''

'''
n=0
for particlenumber in [200]:
    n=+1
    print n
    boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/200)+(4.0/3.0*np.pi))**(1/3.0)
    #complexa=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.5,length,boxsize,particlenumber,0.05)
    #complexa.observables=["number_of_particle"]
    #complexa.show_observable()

    complexb=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,1.0,length,boxsize,particlenumber,0.05)
    complexb.observables=["radial"]
    complexb.show_observable()
'''




for particlenumber1 in [80]:
    lambda_plus1=1
    lambda_c1=0.1
    lambda_minus1=0.1
    length1=2**15
    boxsize1=(((7**3-((4./3.0)*np.pi))*particlenumber1/200)+(4.0/3.0*np.pi))**(1/3.0)
    #complexa=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.5,length,boxsize,particlenumber,0.05)
    #complexa.observables=["number_of_particle"]
    #complexa.show_observable()

    complexb=comlex_anal.Analy_Complex(D,R,lambda_plus1,lambda_minus1,lambda_c1,0.5,length1,boxsize1,particlenumber1,0.05)
    complexb.observables=["number_of_particle"]
    complexb.show_observable()

    # Paramter
    boxsize1=complexb.boxsize*1.0
    S_0=(complexb.particles*1.0)/(boxsize1**3)
    E_0=1.0/(boxsize1**3)
    k2=complexb.k_minus
    k_1=complexb.k_plus
    k_complex=complexb.k_complex
    K_m=(k2+k_complex)/k_1

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



    print k_1, K_m,k2,k_complex
'''
for particlenumber1 in [200]:
    lambda_plus1=1
    lambda_c1=0.1
    lambda_minus1=0.1
    length1=2**15
    boxsize1=(((7**3-((4./3.0)*np.pi))*particlenumber1/200)+(4.0/3.0*np.pi))**(1/3.0)
    #complexa=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.5,length,boxsize,particlenumber,0.05)
    #complexa.observables=["number_of_particle"]
    #complexa.show_observable()

    complexb=comlex_anal.Analy_Complex(D,R,lambda_plus1,lambda_minus1,lambda_c1,0.5,length1,boxsize1,particlenumber1,0.05)
    complexb.observables=["radial"]
    complexb.show_observable()

    # Paramter
    boxsize1=complexb.boxsize*1.0
    S_0=(complexb.particles*1.0)/(boxsize1**3)
    E_0=1.0/(boxsize1**3)
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

    ######


particlenumber=200
boxsize=7.0
D=1.0/6.0
R=1.0
alpha=0.5
lambda_plus=1.0
lambda_c=1.0
lambda_minus=1.0
length=2**11

for particlenumber in [10,50,100]:
    #length=2**15
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

    ######
'''
plt.ylabel('$ \\frac{N}{N_{max}} $',fontsize=10)
plt.xlabel('t in $\\Gamma$',fontsize=10)
plt.legend(loc=4,fontsize='small')
#plt.savefig('../finalreport/data/checking_stability',dpi=300)
plt.show()

