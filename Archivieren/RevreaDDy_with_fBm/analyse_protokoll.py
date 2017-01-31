# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import analysis
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
        :param D: Diffusion coeficient
        :param R: reaction distance
        :param lambda_plus: microscopic rate forward
        :param lambda_minus: microscopic rate backward
        :param lambda_c: microscopic rate complex
        :param alpha: anomalous diffusion exponent
        :param length: trajectory length
        :param boxsize:
        :param particlenumber:
        :param tau: timestep
        :return: plots scaled concentration over time with a with a fit of MM-machnaism (QSSA for )
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

        """
        S. Schnell and C. Mendoza: Closed form solution for time-dependent enzyme
        kinetics, Journal of Theoretical Biology 187, 207 (1997)
        """
        S_t=K_m*lambertw((S_0/K_m)*np.exp((-k_complex*E_0*complexb.time+S_0)/K_m))
        ES_t=((E_0*S_t)/(K_m+S_t))*(1-np.exp(-(K_m+S_t)*k_1*complexb.time))
        E_t=E_0-ES_t
        P_t=S_0-S_t-ES_t
        """
        Plotting
        """
        plt.plot(complexb.time,ES_t/E_0,linestyle=":" ) #complex
        plt.plot(complexb.time,E_t/E_0,linestyle=":")  #enzyme
        plt.plot(complexb.time,S_t/S_0,linestyle=":")  #substrate
        plt.plot(complexb.time,P_t/S_0,linestyle=":")  #product

def show_k1(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau):
        """
        :param D: Diffusion coeficient
        :param R: reaction distance
        :param lambda_plus: microscopic rate forward
        :param lambda_minus: microscopic rate backward
        :param lambda_c: microscopic rate complex
        :param alpha: anomalous diffusion exponent
        :param length: trajectory length
        :param boxsize:
        :param particlenumber:
        :param tau: timestep
        :return:
        """
        complexb=comlex_anal.Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
        complexb.observables=["number_of_particle"]
        complexb.show_observable(show_plot="off")
        boxsize1=complexb.boxsize*1.0
        S_0=(complexb.particles*1.0)/(boxsize1**3)
        E_0=1.0/(boxsize1**3)
        k2=complexb.k_minus
        k_1=complexb.k_plus #*(complexb.time*0.05)**(np.log(complexb.alpha)/np.pi)
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
        complexb.observables=["kplus_neu"]
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
        plt.plot(complexb.time[1:],6*((complexb.time[1:])**complexb.alpha)*complexb.Diffusion,linestyle="--")
'''
plt.figure(figsize=(5,4))

for alpha in [1.5,1.6]:

    #alpha=1.0
    lambda_c=1.0
    lambda_minus=1.0
    #boxsize=15
    show_particle_alpha(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau,2.5)
plt.ylabel('$c/c_{max}$', fontsize=12)
plt.xlabel('t/$\Delta t$', fontsize=12)
plt.legend(loc=3)
plt.fill_between(np.linspace(100,18000),0.8*10**(-5), 10**(-5),color='blue', alpha=0.5)
plt.text(10000, 0.1*10**(-4), 'radial distribution',
        verticalalignment='bottom', horizontalalignment='right', fontsize=15)
#plt.savefig('./Abschlussarbeit/data/normal-scenario-concentrations',dpi=300)
plt.show()


plt.figure(figsize=(5,4))

for lambda_plus in [2.5]:

    alpha=1.0
    lambda_c=1.0
    lambda_minus=0.0
    #boxsize=15
    show_particle_alpha(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau,2.5)
plt.ylabel('$c/c_{max}$', fontsize=12)
plt.xlabel('t/$\Delta t$', fontsize=12)
plt.legend(loc=3)
plt.fill_between(np.linspace(100,18000),0.8*10**(-5), 10**(-5),color='blue', alpha=0.5)
plt.text(10000, 0.1*10**(-4), 'radial distribution',
        verticalalignment='bottom', horizontalalignment='right', fontsize=15)
#plt.savefig('./Abschlussarbeit/data/normal-scenario-concentrations',dpi=300)
plt.show()

for lambda_minus in [0.0]:
    lambda_plus=0.1

    #boxsize=15
    show_particle_alpha(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau,0)
plt.ylabel('$ \\langle c_{\\beta} \\rangle /c_{s}(0);\langle c_{ES}\\rangle /c_{enz}(0)$', fontsize=14)
plt.xlabel('$t/\Delta t$', fontsize=14)
plt.legend(loc=4,fontsize=14)
plt.ylim(10**(-6),10**1)
plt.fill_between(np.linspace(100,18000),0.8*10**(-4), 10**(-4),color='blue', alpha=0.5)
#plt.text(10000, 0.1*10**(-4), 'radial distribution',
 #       verticalalignment='bottom', horizontalalignment='right', fontsize=14)
plt.savefig('./Abschlussarbeit/data/chapman-limit-concentrations1',dpi=300)
plt.show()
for lambda_minus in [1.0]:
    lambda_plus=1.0
    lambda_c=0.1

    #boxsize=15
    show_particle_alpha(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau,0)
for lambda_minus in [0.0]:
    lambda_plus=1.0
    lambda_c=1.0

    #boxsize=15
    show_particle_alpha(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau,0)
plt.ylabel('$ \\langle c_{\\beta} \\rangle /c_{s}(0);\langle c_{ES}\\rangle /c_{enz}(0)$', fontsize=14)
plt.xlabel('$t/\Delta t$', fontsize=14)
plt.legend(loc=4,fontsize=14)
plt.ylim(10**(-6),10**1)
plt.fill_between(np.linspace(100,18000),0.8*10**(-4), 10**(-4),color='blue', alpha=0.5)
#plt.text(10000, 0.1*10**(-4), 'radial distribution',
 #       verticalalignment='bottom', horizontalalignment='right', fontsize=14)
plt.savefig('./Abschlussarbeit/data/highlocaldensityconcentrationversion2',dpi=300)
plt.show()

for lambda_minus in [1.0]:
    lambda_plus=1.0
    lambda_c=1.0

    #boxsize=15
    show_particle_alpha(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau,0)
for lambda_minus in [1.0]:
    lambda_plus=1.0
    alpha=0.5

    #boxsize=15
    show_particle_alpha(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau,1.0)
plt.ylabel('$ \\langle c_{\\beta} \\rangle /c_{s}(0);\langle c_{ES}\\rangle /c_{enz}(0)$', fontsize=14)
plt.xlabel('$t/\Delta t$', fontsize=14)
plt.legend(loc=4,fontsize=14)
plt.ylim(10**(-6),10**1)
plt.fill_between(np.linspace(100,18000),0.8*10**(-4), 10**(-4),color='blue', alpha=0.5)
#plt.text(10000, 0.1*10**(-4), 'radial distribution',
 #       verticalalignment='bottom', horizontalalignment='right', fontsize=14)
plt.savefig('./Abschlussarbeit/data/highlocaldensityconcentrationversionfbmvsbm',dpi=300)
plt.show()

for lambda_minus in [0.0]:
    lambda_plus=1.0
    lambda_c=0.1

    #boxsize=15
    show_particle_alpha(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau,0)
plt.ylabel('$C/C_{max}$', fontsize=14)
plt.xlabel('t', fontsize=14)
plt.legend(loc=2)
plt.fill_between(np.linspace(100/20,18000/20),0.8*10**(-5), 10**(-5),color='blue', alpha=0.5)
plt.text(10000/20, 0.1*10**(-4), 'radial distribution',
        verticalalignment='bottom', horizontalalignment='right', fontsize=15)
#plt.savefig('./Abschlussarbeit/data/chapman-limit-concentrations1',dpi=300)
plt.show()

for lambda_minus in [0.0]:
    lambda_plus=0.1
    lambda_c=1.0
    #boxsize=15
    show_k1_better(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
plt.ylabel('$k_{+}(t)/ \quad \sigma^3/ \Gamma$', fontsize=14)
plt.xlabel('$t/ \Delta t$', fontsize=14)
plt.legend(loc=3,fontsize='small')
#plt.k1_better(np.linspace(100/20,18000/20),0.8*10**(-5), 10**(-5),color='blue', alpha=0.5)

for lambda_minus in [0.0]:
    lambda_plus=1.0
    lambda_c=0.1
    #boxsize=15
    show_k1_better(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
plt.ylabel('$k_{+}(t)/ \quad \sigma^3/ \Gamma$', fontsize=14)
plt.xlabel('$t/ \Delta t$', fontsize=14)
plt.legend(loc=3,fontsize='small')
#plt.savefig('./Abschlussarbeit/data/chapman-limit-concentrations1_k1',dpi=300)

plt.show()
for lambda_plus in [1.0]:
    lambda_minus=1.0
    #lambda_plus=3.0
    lambda_c=0.1
    show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
plt.ylabel('$g(r)$', fontsize=12)
plt.xlabel('r', fontsize=12)
plt.legend(loc=1,fontsize='small')
plt.savefig('./Abschlussarbeit/data/limit-radial-highlocalconcentration',dpi=300)
plt.show()

for lambda_plus in [3.0]:
    lambda_minus=0.0
    #lambda_plus=3.0
    lambda_c=1.0
    show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
plt.ylabel('$g(r)$', fontsize=12)
plt.xlabel('r', fontsize=12)
plt.legend(loc=4,fontsize='small')
plt.savefig('./Abschlussarbeit/data/limit-radial',dpi=300)
plt.show()

for lambda_plus in [3.0]:
    lambda_minus=0.0
    #lambda_plus=3.0
    lambda_c=1.0
    show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
plt.ylabel('$g(r)$', fontsize=12)
plt.xlabel('r', fontsize=12)
plt.legend(loc=4,fontsize='small')
plt.savefig('./Abschlussarbeit/data/limit-radial',dpi=300)
plt.show()

for boxsize in [8.0]:
    lambda_minus=0.0
    lambda_plus=0.1
    show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
plt.ylabel('$g(r)$', fontsize=12)
plt.xlabel('r', fontsize=12)
plt.legend(loc=3,fontsize='small')
plt.savefig('./finalreport/data/chapman-limit-radial',dpi=300)
plt.show()

for alpha in [1.0,0.9,0.8,0.7,0.6,0.5]:
    if alpha==1:
        a=show_particle_alpha(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau,0)
    else:
        show_particle_alpha(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau,a)
plt.ylabel('scaled concentrations', fontsize=10)
plt.xlabel('scaled time', fontsize=10)
plt.legend(loc=3,fontsize='small')
#plt.savefig('./finalreport/data/chapman-limit-radial',dpi=300)
plt.show()

for alpha in [1.0,0.9,0.8,0.7,0.6]:
       show_MSD(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)

plt.ylabel('MSD(t)', fontsize=10)
plt.xlabel('scaled time', fontsize=10)
plt.legend(loc=2,fontsize='small')
#plt.savefig('./finalreport/data/chapman-limit-radial',dpi=300)
plt.show()



for alpha in [1.0,0.9,0.8,0.7,0.6]:
    show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
plt.ylabel('$g(r)$', fontsize=12)
plt.xlabel('r', fontsize=10)
plt.legend(loc=3,fontsize='small')
#plt.savefig('./finalreport/data/chapman-limit-radial_alpha',dpi=300)
plt.show()

plt.figure(figsize=(8,4))
for alpha in [1.2,1.6]:
    show_k1_better(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
#plt.xscale('log')
#plt.yscale('log')
plt.ylabel('$k_{+}(t)/ \quad \sigma^3/ \Gamma$', fontsize=14)
plt.xlabel('$t/ \Delta t$', fontsize=14)
lgd=plt.legend(loc=3,fontsize='small')
plt.savefig('./Abschlussarbeit/data/fractionalalphachangek2neu',bbox_extra_artists=(lgd,), bbox_inches='tight',dpi=300)   # save the figure to file

plt.show()
'''

plt.figure(figsize=(8,4))
for alpha in [1.0]:
    lambda_plus=1.0
    lambda_c=0.1
    lambda_minus=0.0
    show_k1_better(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
#plt.xscale('log')
#plt.yscale('log')
plt.ylabel('$k_{+}(t)/ \quad \sigma^3/ \Gamma$', fontsize=14)
plt.xlabel('$t/ \Delta t$', fontsize=14)
lgd=plt.legend(loc=3,fontsize='small')
#plt.savefig('./Abschlussarbeit/data/fractionalalphachangek2neu',bbox_extra_artists=(lgd,), bbox_inches='tight',dpi=300)   # save the figure to file

plt.show()
'''
lambda_minus=0.0
lambda_plus=3.0
lambda_c=3.0

show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)

lambda_minus=1.0
lambda_plus=1.0
lambda_c=0.1
plt.figure(figsize=(8,4))
show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
plt.ylabel('$g(r)$', fontsize=14)
plt.yscale('log')
plt.xlabel('$r/ \sigma$', fontsize=14)
lgd=plt.legend(loc=1,fontsize='small')
plt.savefig('./Abschlussarbeit/data/limit-radial-highlocalconcentration',bbox_extra_artists=(lgd,), bbox_inches='tight',dpi=300)   # save the figure to file

plt.show()


plt.figure(figsize=(8,4))

lambda_minus=1.0
lambda_plus=1.0
lambda_c=1.0

show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
alpha=0.5
lambda_minus=1.0
lambda_plus=1.0
lambda_c=1.0

show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)

plt.ylabel('$g(r)$', fontsize=14)
plt.xlabel('$r/ \\sigma$', fontsize=14)
plt.yscale('log')
lgd=plt.legend(loc=4,fontsize='small')
plt.savefig('./Abschlussarbeit/data/chapman-limit-radial_alpha_betterfbm',bbox_extra_artists=(lgd,), bbox_inches='tight',dpi=300)   # save the figure to file

plt.show()


lambda_minus=0.0
lambda_c=1.0
lambda_plus=0.1
show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
plt.ylabel('$g(r)$', fontsize=12)
plt.xlabel('r', fontsize=10)
plt.legend(loc=4,fontsize='small')
#plt.savefig('./finalreport/data/chapman-limit-radial_alpha',dpi=300)
plt.show()


for alpha in [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]:
    show_k1_better(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
plt.ylabel('$g(r)$', fontsize=12)
plt.xlabel('r', fontsize=10)
plt.legend(loc=4,fontsize='small')
#plt.savefig('./finalreport/data/chapman-limit-radial_alpha',dpi=300)
plt.show()


plt.figure(figsize=(8,4))
for alpha in [1.0,0.5]:
    lambda_c=0.02
    if alpha==1.0:
        a=show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
    else:
        show_radial(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
plt.ylabel('$g(r)$', fontsize=14)
plt.xlabel('$r/ \sigma$', fontsize=14)
plt.yscale('log')
lgd=plt.legend(loc=1,fontsize=12)
plt.savefig('./Abschlussarbeit/data/chapman-limit-radial_alpha_betterfbmfast',bbox_extra_artists=(lgd,), bbox_inches='tight',dpi=300)   # save the figure to file
plt.show()


for alpha in [0.2]:
       show_MSD(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)

plt.ylabel('MSD(t)', fontsize=10)
plt.xlabel('scaled time', fontsize=10)
plt.legend(loc=2,fontsize='small')
#plt.savefig('./finalreport/data/chapman-limit-radial',dpi=300)
plt.show()


for lambda_plus in [1.0,2.0,3.0]:
    lambda_minus=0.0
    lambda_c=1.0
    alpha=1.0
    length=2**14
    k1_alpha_0=0.0
    show_k1_better(D,R,lambda_plus,lambda_minus,lambda_c,alpha,length,boxsize,particlenumber,tau)
plt.legend(loc=3,fontsize='small')
#plt.savefig('./finalreport/data/chapman-limit-radial',dpi=300)
plt.show()
'''
