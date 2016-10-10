import comlex_anal
import complex_sim
import numpy as np


particlenumber=20
boxsize=8.0
D=1.0/6.0
R=1.0
alpha=1.0
lambda_plus=1.0
lambda_c=1.0
lambda_minus=1.0
length=2**14

for boxsize in [5,6,7,8,9,10,11]:
    for index1 in range(400):
        complex1=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,1.0,length,boxsize,particlenumber,0.05)
        complex1.observables=["all"]
        print index1
        complex1.run()


for lambda_plus in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]:
    #boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/20)+(4.0/3.0*np.pi))**(1/3.0)
    P_in_volume=(1-(4/3*np.pi*R)/(1.0*boxsize**3))**particlenumber
    print P_in_volume
    print boxsize ,"boxsize"
    for index1 in range(400):
        complex=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,1.0,length,boxsize,particlenumber,0.05)
        complex.observables=["all"]
        print R
        print index1
        complex.run()

for lambda_minus in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.00]:
    #boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/20)+(4.0/3.0*np.pi))**(1/3.0)
    P_in_volume=(1-(4/3*np.pi*R)/(1.0*boxsize**3))**particlenumber
    print P_in_volume
    print boxsize ,"boxsize"
    for index1 in range(400):
        complex=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,1.0,length,boxsize,particlenumber,0.05)
        complex.observables=["all"]
        print R
        print index1
        complex.run()

for lambda_c in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.00]:
    #boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/20)+(4.0/3.0*np.pi))**(1/3.0)
    P_in_volume=(1-(4/3*np.pi*R)/(1.0*boxsize**3))**particlenumber
    print P_in_volume
    print boxsize ,"boxsize"
    for index1 in range(400):
        complex=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,1.0,length,boxsize,particlenumber,0.05)
        complex.observables=["all"]
        print R
        print index1
        complex.run()
for boxsize in [5,6,7,8,9,10,11]:
    #boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/20)+(4.0/3.0*np.pi))**(1/3.0)
    P_in_volume=(1-(4/3*np.pi*R)/(1.0*boxsize**3))**particlenumber
    print P_in_volume
    print boxsize ,"boxsize"

    for index1 in range(400):

        complex=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.5,length,boxsize,particlenumber,0.05)

        complex.observables=["all"]
        print R
        print index1
        complex.run()

for lambda_plus in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.00]:
    #boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/20)+(4.0/3.0*np.pi))**(1/3.0)
    P_in_volume=(1-(4/3*np.pi*R)/(1.0*boxsize**3))**particlenumber
    print P_in_volume
    print boxsize ,"boxsize"
    for index1 in range(400):
        complex=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.5,length,boxsize,particlenumber,0.05)
        complex.observables=["all"]
        print R
        print index1
        complex.run()

for lambda_minus in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.00]:
    #boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/20)+(4.0/3.0*np.pi))**(1/3.0)
    P_in_volume=(1-(4/3*np.pi*R)/(1.0*boxsize**3))**particlenumber
    print P_in_volume
    print boxsize ,"boxsize"
    for index1 in range(400):
        complex=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.5,length,boxsize,particlenumber,0.05)
        complex.observables=["all"]
        print R
        print index1
        complex.run()

for lambda_c in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.00]:
    #boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/20)+(4.0/3.0*np.pi))**(1/3.0)
    P_in_volume=(1-(4/3*np.pi*R)/(1.0*boxsize**3))**particlenumber
    print P_in_volume
    print boxsize ,"boxsize"
    for index1 in range(400):
        complex=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.5,length,boxsize,particlenumber,0.05)
        complex.observables=["all"]
        print R
        print index1
        complex.run()

for alpha in [0.5,0.6,0.7,0.8,0.9,1.00]:
    #boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/20)+(4.0/3.0*np.pi))**(1/3.0)
    P_in_volume=(1-(4/3*np.pi*R)/(1.0*boxsize**3))**particlenumber
    print P_in_volume
    print boxsize ,"boxsize"
    for index1 in range(400):
        complex=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,0.05)
        complex.observables=["all"]
        print R
        print index1
        complex.run()



    for index2 in range(400):
        complex1=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,1.0,length,boxsize,particlenumber,0.05)
        complex1.observables=["all"]
        #print index2
        #complex1.run()



'''
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
    '''


'''
plt.ylabel('radial distribution function',fontsize=10)
plt.xlabel('r',fontsize=10)
plt.legend(loc=4,fontsize='small')
#plt.savefig('../finalreport/data/checking_stability',dpi=300)
plt.show()

'''
