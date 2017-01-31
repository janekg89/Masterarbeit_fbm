import logging
import os
import random
import string
import time
#from __future__ import print_function
import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
import simulation_model
from sklearn import datasets, linear_model
from scipy.optimize import curve_fit
from scipy.special import lambertw
from scipy import stats

reload(logging)
logging.basicConfig(
    format='%(asctime)s %(levelname)s: %(message)s',
    datefmt='%Y/%m/%d %I:%M:%S',
    level=logging.DEBUG
    )
class Analy_Complex(simulation_model.Sim_Complex):
    def func(self,params, X, Y):
                    # extract current values of fit parameters from input array
                    K_m = params[0]
                    k2=params[1]
                    # compute chi-square
                    chi2 = 0.0
                    S_0=(self.particles*1.0)/(1.0*self.boxsize**3)
                    E_0=1.0/(1.0*self.boxsize**3)
                            #k2=self.micro_reactionrate_complex
                    for n in range(len(X)):
                        x = X[n]
                        y=K_m*lambertw((S_0/K_m)*np.exp((-k2*E_0*x+S_0)/K_m))
                        chi2 = chi2 + (Y[n] - y)**2
                    return chi2
    def func_es(self,params, X, Y):
                    # extract current values of fit parameters from input array
                    K_m = params[0]
                    k2=params[1]
                    # compute chi-square
                    chi2 = 0.0
                    S_0=(self.particles*1.0)/(1.0*self.boxsize**3)
                    E_0=1.0/(1.0*self.boxsize**3)
                    for n in range(len(X)):
                        x = X[n]
                        y_old=K_m*lambertw((S_0/K_m)*np.exp((-k2*E_0*x+S_0)/K_m))
                        y=(((E_0*y_old)/(K_m+y_old))*(1-np.exp(-(K_m+y_old)*((self.micro_reactionrate_backward+self.micro_reactionrate_complex)/K_m)*x)))
                        chi2 = chi2 + (Y[n] - y)**2
                    return chi2
    def func_p(self,params, X, Y):
                    K_m = params[0]

                    chi2 = 0.0
                    S_0=(self.particles*1.0)/(1.0*self.boxsize**3)
                    E_0=1.0/(1.0*self.boxsize**3)
                    k2=self.micro_reactionrate_complex
                    for n in range(len(X)):
                        x = X[n]
                        y_old=K_m*lambertw((S_0/K_m)*np.exp((-k2*E_0*x+S_0)/K_m))
                        y=(((E_0*y_old)/(K_m+y_old))*(1-np.exp(-(K_m+y_old)*((self.micro_reactionrate_backward+self.micro_reactionrate_complex)/K_m)*x)))
                        y_neu=S_0-y_old-y
                        chi2 = chi2 + (Y[n] - y_neu)**2


                    return chi2
    def show_observable(self,show_plot="on",k1_alpha_0=0):
        self.k1_alpha_0=k1_alpha_0
        for ii in self.observables:

            if ii == "number_of_particle":
                particle_number_c=[]
                particle_number_s=[]
                particle_number_p=[]
                for file in os.listdir(self.path):
                    if  file.endswith("number_c.h5"):
                        f = h5.File(self.path+file, "r")
                        t = f['times'][:]
                        fc = f['particleNumbers'][:]
                        particle_number_c.append(fc)
                        f.close()
                    elif  file.endswith("number_s.h5"):
                        f = h5.File(self.path+file, "r")
                        t = f['times'][:]
                        fc = f['particleNumbers'][:]
                        particle_number_s.append(fc)
                        f.close()
                    elif  file.endswith("number_p.h5"):
                        f = h5.File(self.path+file, "r")
                        t = f['times'][:]
                        fc = f['particleNumbers'][:]
                        particle_number_p.append(fc)
                        f.close()


                self.time=t
                particle_number_c=np.array(particle_number_c)
                #print particle_number_c.shape
                #print t[0:10]
                #print self.timestep
                particle_number_c_mean=particle_number_c.mean(axis=0)
                linestyle='-'
                marker=""
                markersize=0.5
                every=1
                color="purple"
                if self.alpha==0.5:
                    color="b"
                    linestyle='--'
                    marker=""
                    #markersize=1.0

                plt.plot(t, particle_number_c_mean,marker=marker,linestyle=linestyle,color=color,markevery=every ,label="complex") #of complex

                particle_number_s=np.array(particle_number_s)
                particle_number_s_mean=particle_number_s.mean(axis=0)
                if show_plot is "on":
                    plt.plot(t, particle_number_s_mean/(self.particles*1.0),linestyle=linestyle,color="r" ,marker=marker,markevery=every, label="substrate") #of substrate

                particle_number_p=np.array(particle_number_p)
                particle_number_p_mean=particle_number_p.mean(axis=0)


                def s_t_modell(t,kplus,kminus,kcomplex):
                    S_0=(self.particles*1.0)/(self.boxsize**3)
                    E_0=1.0/(self.boxsize**3)
                    kcomplex=self.micro_reactionrate_complex
                    kminus=self.micro_reactionrate_backward
                    #t= np.linspace(0,self.length*self.timestep,10000)
                    S_t=((kminus+kcomplex)/kplus)*lambertw((S_0/((kminus+kcomplex)/kplus))*np.exp((-kcomplex*E_0*t+S_0)/((kminus+kcomplex)/kplus)))
                    #S_t=((kminus+kcomplex)/(kplus*(t*0.05)**(np.log(self.alpha)/np.pi)))*lambertw((S_0/((kminus+kcomplex)/(kplus*(t*0.05)**(np.log(self.alpha)/np.pi))))*np.exp((-kcomplex*E_0*t+S_0)/((kminus+kcomplex)/(kplus*(t*0.05)**(np.log(self.alpha)/np.pi)))))
                    return S_t.real
                def es_t_modell(t,kplus,kminus,kcomplex):
                    #kminus=1
                    kcomplex=self.micro_reactionrate_complex
                    kminus=self.micro_reactionrate_backward
                    S_0=(self.particles*1.0)/(self.boxsize**3)
                    E_0=1.0/(self.boxsize**3)
                    #t= np.linspace(0,self.length*self.timestep,10000)
                    S_t=((kminus+kcomplex)/kplus)*lambertw((S_0/((kminus+kcomplex)/kplus))*np.exp((-kminus*E_0*t+S_0)/((kminus+kcomplex)/kplus)))
                    ES_t=((E_0*S_t.real)/(((kminus+kcomplex)/kplus)+S_t.real))*(1-np.exp(-(((kminus+kcomplex)/kplus)+S_t.real)*kplus*t))

                    return ES_t
                def p_t_modell(t,kplus,kminus,kcomplex):
                    #kminus=1
                    kcomplex=self.micro_reactionrate_complex
                    kminus=self.micro_reactionrate_backward
                    S_0=(self.particles*1.0)/(self.boxsize**3)
                    E_0=1.0/(self.boxsize**3)
                    #t= np.linspace(0,self.length*self.timestep,10000)
                    S_t=((kminus+kcomplex)/kplus)*lambertw((S_0/((kminus+kcomplex)/kplus))*np.exp((-kminus*E_0*t+S_0)/((kminus+kcomplex)/kplus)))
                    ES_t=((E_0*S_t.real)/(((kminus+kcomplex)/kplus)+S_t.real))*(1-np.exp(-(((kminus+kcomplex)/kplus)+S_t.real)*kplus*t))
                    P_t=S_0-S_t.real-ES_t
                    return P_t

                best_vals2, covar2= curve_fit(p_t_modell,t,(particle_number_p_mean*1.0)/(1.0*self.boxsize**3))
                best_vals1, covar1= curve_fit(s_t_modell,t,(particle_number_s_mean*1.0)/(1.0*self.boxsize**3))
                best_vals, covar= curve_fit(es_t_modell,t,(particle_number_c_mean*1.0)/(1.0*self.boxsize**3))


                print best_vals[0],best_vals1[0],best_vals2[0]

                self.k_all=[best_vals[0],best_vals1[0],best_vals2[0]]
                self.var_all=[covar[0],covar1[0],covar2[0]]

                self.k_plus=best_vals1[0]
                self.k_minus=self.micro_reactionrate_backward
                self.k_complex=self.micro_reactionrate_complex
                plt.plot(t, particle_number_p_mean/(self.particles*1.0),linestyle=linestyle,color="g"  ,marker=marker,markevery=every,label=" product " ) # of product


            if ii == "Michaelis_fit_alpha":
                particle_number_c=[]
                particle_number_s=[]
                particle_number_p=[]

                reaction1_f=[]
                for file in os.listdir(self.path+"number_C/"):
                    f = h5.File(self.path+"number_C/"+file, "r")
                    particle_number_c.append(f['particleNumbers'][:])
                    f.close()
                for file in os.listdir(self.path+"number_S/"):
                        f = h5.File(self.path+"number_S/"+file, "r")
                        particle_number_s.append(f['particleNumbers'][:])
                        f.close()
                for file in os.listdir(self.path+"number_P/"):
                        f = h5.File(self.path+"number_P/"+file, "r")
                        t = f['times'][:]
                        particle_number_p.append(f['particleNumbers'][:])
                        f.close()
                for file in os.listdir(self.path+"reaction1/"):
                        f = h5.File(self.path+"reaction1/"+file, "r")
                        reaction1_f.append(f['forwardCounter'][:])
                        t = f['times'][:]
                        f.close()



                self.time=t

                particle_number_s=np.array(particle_number_s)
                print particle_number_s.shape, "das ist die Anzahl der Trajektorien mal die Laenge der Trajektorien"
                particle_number_s_mean=particle_number_s.mean(axis=0)
                particle_number_c=np.array(particle_number_c)
                particle_number_c_mean=particle_number_c.mean(axis=0)
                particle_number_p=np.array(particle_number_p)
                particle_number_p_mean=particle_number_p.mean(axis=0)

                reaction1_f=np.array(reaction1_f)
                reaction1_f_mean=reaction1_f.mean(axis=0)

                k1=(np.diff(reaction1_f_mean)*self.boxsize**3)/((particle_number_s_mean[1:])*(1-particle_number_c_mean[1:]))


                linestyle='-'
                marker="+"
                markersize=0.5
                every=1
                color="purple"


                    #markersize=1.0


                def s_t_modell(t,kplus,kminus,kcomplex):
                    S_0=(self.particles*1.0)/(self.boxsize**3)
                    E_0=1.0/(self.boxsize**3)
                    kcomplex=self.micro_reactionrate_complex
                    kminus=self.micro_reactionrate_backward
                    S_t=((kminus+kcomplex)/kplus)*lambertw((S_0/((kminus+kcomplex)/kplus))*np.exp((-kcomplex*E_0*t+S_0)/((kminus+kcomplex)/kplus)))
                    return S_t.real
                def es_t_modell(t,kplus,kminus,kcomplex):
                    kcomplex=self.micro_reactionrate_complex
                    kminus=self.micro_reactionrate_backward
                    S_0=(self.particles*1.0)/(self.boxsize**3)
                    E_0=1.0/(self.boxsize**3)
                    S_t=((kminus+kcomplex)/kplus)*lambertw((S_0/((kminus+kcomplex)/kplus))*np.exp((-kcomplex*E_0*t+S_0)/((kminus+kcomplex)/kplus)))
                    ES_t=((E_0*S_t.real)/(((kminus+kcomplex)/kplus)+S_t.real))*(1-np.exp(-(((kminus+kcomplex)/kplus)+S_t.real)*kplus*t))
                    return ES_t

                def p_t_modell(t,kplus,kminus,kcomplex):
                    kcomplex=self.micro_reactionrate_complex
                    kminus=self.micro_reactionrate_backward
                    S_0=(self.particles*1.0)/(self.boxsize**3)
                    E_0=1.0/(self.boxsize**3)
                    S_t=((kminus+kcomplex)/kplus)*lambertw((S_0/((kminus+kcomplex)/kplus))*np.exp((-kcomplex*E_0*t+S_0)/((kminus+kcomplex)/kplus)))
                    ES_t=((E_0*S_t.real)/(((kminus+kcomplex)/kplus)+S_t.real))*(1-np.exp(-(((kminus+kcomplex)/kplus)+S_t.real)*kplus*t))
                    P_t=S_0-S_t.real-ES_t
                    return P_t



                best_vals2, covar2= curve_fit(p_t_modell,t,(particle_number_p_mean*1.0)/(1.0*self.boxsize**3))
                best_vals1, covar1= curve_fit(s_t_modell,t,(particle_number_s_mean*1.0)/(1.0*self.boxsize**3))
                best_vals, covar= curve_fit(es_t_modell,t,(particle_number_c_mean*1.0)/(1.0*self.boxsize**3))


                self.k_plus=best_vals2[0]
                print best_vals[0],best_vals1[0],best_vals2[0]

                def makemean(k1,step):
                    i=0
                    index_old=0
                    t_neu=[]
                    t_increment=((step+1)/2)
                    k1_neu=[]
                    while i<self.length/step:
                        index=index_old+step
                        k1_neu.append(k1[index_old:index].mean())
                        t_neu.append(t_increment)
                        t_increment=t_increment+step
                        index_old=index
                        i=i+1
                    return np.array(t_neu), np.array(k1_neu)

                t_index,k1neu=makemean(k1,3)
                tneu=t[t_index]
                if self.alpha==1.0:
                    self.k1_alpha_0=self.k_plus
                print self.k1_alpha_0



                assert self.k1_alpha_0 != 0.0,"put k1 of alpha = 0 as an argument of show_observable"


                def k1_fit(t,h):

                    return self.k1_alpha_0*t**(-h)


                h_best, covar_h= curve_fit(k1_fit,tneu,k1neu*20)
                t_index,k1neu=makemean(k1,151)
                tneu=t[t_index]



                print h_best,"das ist h"
                print self.k1_alpha_0,"k_0"
                print self.alpha, "das ist alpha"




                def s_t_modell_alpha(t,h,kminus,kcomplex):
                    kplus=self.k1_alpha_0*t**(-h)
                    S_0=(self.particles*1.0)/(self.boxsize**3)
                    E_0=1.0/(self.boxsize**3)
                    kcomplex=self.micro_reactionrate_complex
                    kminus=self.micro_reactionrate_backward
                    S_t=((kminus+kcomplex)/kplus)*lambertw((S_0/((kminus+kcomplex)/kplus))*np.exp((-kcomplex*E_0*t+S_0)/((kminus+kcomplex)/kplus)))
                    return S_t.real

                def p_t_modell_alpha(t,h,kminus,kcomplex):

                    kplus=self.k1_alpha_0*t**(-h)
                    kcomplex=self.micro_reactionrate_complex
                    kminus=self.micro_reactionrate_backward
                    S_0=(self.particles*1.0)/(self.boxsize**3)
                    E_0=1.0/(self.boxsize**3)
                    S_t=((kminus+kcomplex)/kplus)*lambertw((S_0/((kminus+kcomplex)/kplus))*np.exp((-kcomplex*E_0*t+S_0)/((kminus+kcomplex)/kplus)))
                    ES_t=((E_0*S_t.real)/(((kminus+kcomplex)/kplus)+S_t.real))*(1-np.exp(-(((kminus+kcomplex)/kplus)+S_t.real)*kplus*t))
                    P_t=S_0-S_t.real-ES_t
                    return P_t
                def es_t_modell_alpha(t,h,kminus,kcomplex):
                    kplus=self.k1_alpha_0*t**(-h)
                    kcomplex=self.micro_reactionrate_complex
                    kminus=self.micro_reactionrate_backward
                    S_0=(self.particles*1.0)/(self.boxsize**3)
                    E_0=1.0/(self.boxsize**3)
                    S_t=((kminus+kcomplex)/kplus)*lambertw((S_0/((kminus+kcomplex)/kplus))*np.exp((-kcomplex*E_0*t+S_0)/((kminus+kcomplex)/kplus)))
                    ES_t=((E_0*S_t.real)/(((kminus+kcomplex)/kplus)+S_t.real))*(1-np.exp(-(((kminus+kcomplex)/kplus)+S_t.real)*kplus*t))
                    return ES_t



                print self.k1_alpha_0,"k_0"
                print self.alpha, "das ist alpha"





                S_0=(self.particles*1.0)/(1*self.boxsize**3)
                E_0=1.0/(1*self.boxsize**3)

                marker="+"
                '''
                if self.micro_reactionrate_complex==0.1:
                    plt.loglog(t[1:]/self.timestep, particle_number_c_mean[1:],color=color,marker=marker,linestyle="",markevery=np.logspace(0,np.log10(len(range(2**14))),60).astype('int')[2:] ,label="complex ") #of complex
                    plt.plot(t[1:]/self.timestep, particle_number_s_mean[1:]/(self.particles*1.0),linestyle=linestyle,color="r" , label="substrate ") #of substrate
                    plt.plot(t[1:]/self.timestep, particle_number_p_mean[1:]/(self.particles*1.0),linestyle=linestyle,color="g"  ,label=" product " ) # of product
                '''
                tneu=np.array(range(2**16-1))

                if show_plot is "on":
                    if self.alpha==1.0:
                        plt.loglog(t[np.logspace(0,np.log10(len(range(2**14))),60).astype('int')[2:]]/self.timestep, particle_number_c_mean[np.logspace(0,np.log10(len(range(2**14))),60).astype('int')[2:]],color=color ,alpha=0.3) #of complex
                        plt.plot(t[1:]/self.timestep, particle_number_s_mean[1:]/(self.particles*1.0),linestyle=linestyle,color="r" ,alpha=0.3) #of substrate
                        plt.plot(t[1:]/self.timestep, particle_number_p_mean[1:]/(self.particles*1.0),linestyle=linestyle,color="g" ,alpha=0.3  ) # of product
                    else:
                        plt.loglog(t[np.logspace(0,np.log10(len(range(2**14))),60).astype('int')[2:]]/self.timestep, particle_number_c_mean[np.logspace(0,np.log10(len(range(2**14))),60).astype('int')[2:]],color=color ,label="$\\langle c_{ES} (t) \\rangle$ ") #of complex
                        plt.plot(t[1:]/self.timestep, particle_number_s_mean[1:]/(self.particles*1.0),linestyle=linestyle,color="r" , label="$ \\langle c_{S} (t)\\rangle$ ") #of substrate
                        plt.plot(t[1:]/self.timestep, particle_number_p_mean[1:]/(self.particles*1.0),linestyle=linestyle,color="g"  ,label=" $ \\langle c_{P} (t)\\rangle$" ) # of product

                    #plt.plot(tneu/self.timestep,p_t_modell_alpha(tneu,self.k1_alpha_0,1.0,1.0)/S_0,color=color )
                    #plt.plot(tneu/self.timestep,s_t_modell_alpha(tneu,k1_fit(tneu,h_best),1.0,1.0)/S_0,color="r" )
                    #plt.plot(tneu/self.timestep,es_t_modell_alpha(tneu,k1_fit(t[t_index],h_best),1.0,1.0)*(1*self.boxsize**3),color="g" )

                    #if self.alpha==0.5:
                        #plt.scatter(tneu[np.logspace(0,np.log10(len(tneu)),40).astype('int')[2:]]/self.timestep,p_t_modell(tneu, self.k1_alpha_0,self.micro_reactionrate_backward,self.micro_reactionrate_complex)[np.logspace(0,np.log10(len(tneu)),40).astype('int')[2:]]/S_0,alpha=0.3,color="g",facecolors='none' )
                        #plt.scatter(tneu[np.logspace(0,np.log10(len(tneu)),40).astype('int')[2:]]/self.timestep,s_t_modell(tneu,self.k1_alpha_0,self.micro_reactionrate_backward,self.micro_reactionrate_complex)[np.logspace(0,np.log10(len(tneu)),40).astype('int')[2:]]/S_0,alpha=0.3,color="r",facecolors='none')
                        #plt.scatter(tneu[np.logspace(0,np.log10(len(tneu)),40).astype('int')[2:]]/self.timestep,es_t_modell(tneu,self.k1_alpha_0,self.micro_reactionrate_backward,self.micro_reactionrate_complex)[np.logspace(0,np.log10(len(tneu)),40).astype('int')[2:]]*self.boxsize**3,alpha=0.3,color=color,facecolors='none')


                    if self.alpha==1.0:
                        plt.scatter(tneu[np.logspace(0,np.log10(len(tneu)),40).astype('int')[2:]]/self.timestep,p_t_modell(tneu, self.k1_alpha_0,self.micro_reactionrate_backward,self.micro_reactionrate_complex)[np.logspace(0,np.log10(len(tneu)),40).astype('int')[2:]]/S_0,alpha=0.3,color="g",facecolors='none' )
                        plt.scatter(tneu[np.logspace(0,np.log10(len(tneu)),40).astype('int')[2:]]/self.timestep,s_t_modell(tneu,self.k1_alpha_0,self.micro_reactionrate_backward,self.micro_reactionrate_complex)[np.logspace(0,np.log10(len(tneu)),40).astype('int')[2:]]/S_0,alpha=0.3,color="r",facecolors='none')
                        plt.scatter(tneu[np.logspace(0,np.log10(len(tneu)),40).astype('int')[2:]]/self.timestep,es_t_modell(tneu,self.k1_alpha_0,self.micro_reactionrate_backward,self.micro_reactionrate_complex)[np.logspace(0,np.log10(len(tneu)),40).astype('int')[2:]]*self.boxsize**3,alpha=0.3,color=color,facecolors='none')
                """
                name='particle'+'.txt'
                with open(self.path+name, 'w') as particlenumber:
                        particlenumber.write("The Amount of mean particles at time t")
                        particlenumber.write("c;s;e")
                        for i in range(len(particle_number_c_mean)):
                            particlenumber.write("%d;%d;%d,%d" %(particle_number_c_mean[i],particle_number_s_mean[i]/(self.particles*1.0),particle_number_p_mean[i]/(self.particles*1.0),t[i]))

                name='MMmodel'+'.txt'
                with open(self.path+name, 'w') as particlenumber:
                        particlenumber.write("The Amount of mean particles at time t")
                        particlenumber.write("c;s;e")
                        for i in range(len(es_t_modell(tneu,self.k1_alpha_0,1.0,1.0)*(1*self.boxsize**3))):
                            particlenumber.write("%d;%d;%d;%d" %(es_t_modell(tneu,self.k1_alpha_0,1.0,1.0)[i]*(1*self.boxsize**3),s_t_modell(tneu,self.k1_alpha_0,1.0,1.0)[i]/S_0,p_t_modell(tneu,self.k1_alpha_0,1.0,1.0)[i]/S_0,tneu[i]))


                """


                name='mean_values_particle'+'.h5'
                with h5.File(self.path+name, 'w') as hf:
                    hf.create_dataset('norm_number_c', data=particle_number_c_mean)
                    hf.create_dataset('norm_number_s', data=particle_number_s_mean/(self.particles*1.0))
                    hf.create_dataset('norm_number_p', data=particle_number_p_mean/(self.particles*1.0))
                    hf.create_dataset('t', data=t)






            elif ii== "radial":
                radial_s=[]
                radial_p=[]
                radial_sp=[]
                linestyle='-'
                if self.alpha==0.5:
                    linestyle='--'
                for file in os.listdir(self.path+"radialS/"):
                    f = h5.File(self.path+"radialS/"+file, "r")
                    radial_s.append(f['bins'][:])
                    f.close()


                for file in os.listdir(self.path+"radialP/"):
                    f = h5.File(self.path+"radialP/"+file, "r")
                    bincenters = f['binCenters'][:]
                    radial_p.append(f['bins'][:])
                    f.close()

                for file in os.listdir(self.path+"radialSP/"):
                    f = h5.File(self.path+"radialSP/"+file, "r")
                    bincenters = f['binCenters'][:]
                    radial_sp.append(f['bins'][:])

                    f.close()


                radial_s=np.array(radial_s)
                radial_p=np.array(radial_p)
                radial_sp=np.array(radial_sp)
                if self.micro_reactionrate_complex==1.0 and self.micro_reactionrate_backward==0.0 :
                    plt.errorbar(bincenters,(1./1.05)*radial_s.mean(axis=0)*self.boxsize**3/len(self.radial_time_range),stats.sem(radial_s,axis=0)*self.boxsize**3/len(self.radial_time_range),linestyle="",label="$\\lambda_c=%.2f$ $\\lambda_{+}=%.2f$" %(self.micro_reactionrate_complex, self.micro_reactionrate_forward))

                if self.micro_reactionrate_complex==3.0 and self.micro_reactionrate_backward==0.0 :
                    plt.errorbar(bincenters,(1./1.25)*radial_s.mean(axis=0)*self.boxsize**3/len(self.radial_time_range),stats.sem(radial_s,axis=0)*self.boxsize**3/len(self.radial_time_range),linestyle="",label="$\\lambda_c=%.2f$ $\\lambda_{+}=%.2f$" %(self.micro_reactionrate_complex, self.micro_reactionrate_forward))
                if self.micro_reactionrate_complex==0.1 and self.micro_reactionrate_backward==1.0 :
                    plt.errorbar(bincenters,(1./1.05)*radial_s.mean(axis=0)*self.boxsize**3/len(self.radial_time_range),stats.sem(radial_s,axis=0)*self.boxsize**3/len(self.radial_time_range),linestyle="",label="$\\lambda_c=%.2f$ $\\lambda_{+}=%.2f$" %(self.micro_reactionrate_complex, self.micro_reactionrate_forward))
                #plt.
                elif self.alpha==0.5:
                    plt.errorbar(bincenters,(1./1.07)*radial_s.mean(axis=0)*self.boxsize**3/len(self.radial_time_range),stats.sem(radial_s,axis=0)*self.boxsize**3/len(self.radial_time_range),linestyle="",label="$\\alpha=%.2f$ " %(self.alpha))
                elif self.alpha==1.0:
                    plt.errorbar(bincenters,(1./1.04)*radial_s.mean(axis=0)*self.boxsize**3/len(self.radial_time_range),stats.sem(radial_s,axis=0)*self.boxsize**3/len(self.radial_time_range),linestyle="",color='g',label="$\\alpha=%.2f$ " %(self.alpha))


                def radial_erban(bincenters,radial_dist):
                    c_inf=radial_dist[-1]
                    self.Diffusion=self.Diffusion
                    a_2=c_inf*(np.sqrt(self.Diffusion/self.micro_reactionrate_forward)*np.tanh(self.reactiondistance*np.sqrt(self.micro_reactionrate_forward/self.Diffusion))-self.reactiondistance)
                    return c_inf+a_2/bincenters
                def index_cut():
                    for i in range(len(bincenters)):
                        if bincenters[i] >= 1.:
                            return i
                def radial_erban_inside(bincenters,radial_dist):
                    c_inf=radial_dist[-1]
                    a_3=c_inf*np.sqrt(self.Diffusion/self.micro_reactionrate_forward)*(2*np.cosh(self.reactiondistance*np.sqrt(self.micro_reactionrate_forward/self.Diffusion)))**(-1)
                    return ((2*a_3)/bincenters)*np.sinh(bincenters*np.sqrt(self.micro_reactionrate_forward/self.Diffusion))

                #bincentersnew=bincenters[:index_cut()]
                labelweg="Erban-Chapman"
                if self.micro_reactionrate_complex==3.0:
                 labelweg=""



                def fit_radial_ouside(a1,a2,bincenters):
                    return a1 #todo hier von erban chapman



            elif ii== "radial_lambdaplus01_lambdaminus0":
                radial_s=[]
                radial_p=[]
                radial_sp=[]
                linestyle='-'
                if self.alpha==0.5:
                    linestyle='--'
                for file in os.listdir(self.path+"radialS/"):
                    f = h5.File(self.path+"radialS/"+file, "r")
                    radial_s.append(f['bins'][:])
                    f.close()


                for file in os.listdir(self.path+"radialP/"):
                    f = h5.File(self.path+"radialP/"+file, "r")
                    bincenters = f['binCenters'][:]
                    radial_p.append(f['bins'][:])
                    f.close()

                for file in os.listdir(self.path+"radialSP/"):
                    f = h5.File(self.path+"radialSP/"+file, "r")
                    bincenters = f['binCenters'][:]
                    radial_sp.append(f['bins'][:])

                    f.close()


                radial_s=np.array(radial_s)
                radial_p=np.array(radial_p)
                radial_sp=np.array(radial_sp)
                plt.errorbar(bincenters,radial_s.mean(axis=0)*self.boxsize**3/len(self.radial_time_range),stats.sem(radial_s,axis=0)*self.boxsize**3/len(self.radial_time_range),label="$\\lambda_c=%.2f$ $\\lambda_{+}=%.2f$" %(self.micro_reactionrate_complex, self.micro_reactionrate_forward))
                #plt.errorbar(bincenters,radial_p.mean(axis=0)*self.boxsize**3/2,radial_p.std(axis=0)*self.boxsize**3)
                #plt.errorbar(bincenters,radial_sp.mean(axis=0)*self.boxsize**3/len(self.radial_time_range),radial_sp.std(axis=0)*self.boxsize**3)


                def radial_erban(bincenters,radial_dist):
                    c_inf=radial_dist[-1]
                    self.Diffusion=self.Diffusion
                    a_2=c_inf*(np.sqrt(self.Diffusion/self.micro_reactionrate_forward)*np.tanh(self.reactiondistance*np.sqrt(self.micro_reactionrate_forward/self.Diffusion))-self.reactiondistance)
                    return c_inf+a_2/bincenters
                def index_cut():
                    for i in range(len(bincenters)):
                        if bincenters[i] >= 1.:
                            return i
                def radial_erban_inside(bincenters,radial_dist):
                    c_inf=radial_dist[-1]
                    a_3=c_inf*np.sqrt(self.Diffusion/self.micro_reactionrate_forward)*(2*np.cosh(self.reactiondistance*np.sqrt(self.micro_reactionrate_forward/self.Diffusion)))**(-1)
                    return ((2*a_3)/bincenters)*np.sinh(bincenters*np.sqrt(self.micro_reactionrate_forward/self.Diffusion))

                #bincentersnew=bincenters[:index_cut()]

                if self.alpha  == 1.0:
                 plt.plot(bincenters[index_cut():],radial_erban(bincenters,1.15*radial_s.mean(axis=0)*self.boxsize**3/len(self.radial_time_range))[index_cut():],linestyle="--", label="Erban-Chapman")
                 plt.plot(bincenters[:index_cut()+1],radial_erban_inside(bincenters,1.15*radial_s.mean(axis=0)*self.boxsize**3/len(self.radial_time_range))[:index_cut()+1],linestyle="--")
                else:
                 plt.plot(bincenters[index_cut():],radial_erban(bincenters,1.15*radial_s.mean(axis=0)*self.boxsize**3/len(self.radial_time_range))[index_cut():],linestyle="--")
                 plt.plot(bincenters[:index_cut()+1],radial_erban_inside(bincenters,1.15*radial_s.mean(axis=0)*self.boxsize**3/len(self.radial_time_range))[:index_cut()+1],linestyle="--")







                def fit_radial_ouside(a1,a2,bincenters):
                    return a1 #todo hier von erban chapman



            elif ii=="MSD":
                msd_s=[]

                for file in os.listdir(self.path+"msd/"):
                        f = h5.File(self.path+"msd/"+file, "r")
                        msd = f['meanSquaredDisplacements'][:]
                        t = f['times'][:]
                        if len(msd)-1 == self.length:
                            msd_s.append(np.array(msd))

                        f.close()
                msd_s=np.array(msd_s)
                msd_s_mean=msd_s.mean(axis=0)
                print len(msd_s_mean[1:])
                plt.plot(np.array(range(len(msd_s_mean[1:])))/self.timestep,msd_s_mean[1:],label="alpha = %.2f" %self.alpha)

                self.time=t


            elif ii=="reaction":
                reaction1_f=[]
                reaction1_b=[]
                reaction2_b=[]
                S_particle=[]
                ES_particle=[]

                for file in os.listdir(self.path):
                    if  file.endswith("reac_num_S_and_E_to_C.h5"):
                        numberfile=file[:-24]+"number_s.h5"

                        fn = h5.File(self.path+numberfile, "r")
                        particles_s = fn['particleNumbers'][:]

                        numberfile1=file[:-24]+"number_c.h5"

                        fc = h5.File(self.path+numberfile1, "r")
                        particles_c = fc['particleNumbers'][:]

                        S_particle.append(particles_s)
                        ES_particle.append(particles_c)

                        #print np.min(S_concentration)


                        f = h5.File(self.path+file, "r")
                        reaction_f = f['forwardCounter'][:]
                        reaction_b=f['backwardCounter'][:]
                        reaction1_b.append(np.float64(reaction_b))





                        reaction1_f.append(reaction_f)
                        t = f['times'][:]
                        f.close()
                        fn.close()
                        fc.close()


                    elif  file.endswith('reac_num_S_and_P_to_C.h5'):
                        f = h5.File(self.path+file, "r")
                        reaction_b=f['backwardCounter'][:]
                        t = f['times'][:]
                        reaction2_b.append(np.float64(reaction_b))
                        f.close()


                reaction2_b=np.array(reaction2_b)
                reaction2_b_mean=reaction2_b.mean(axis=0)
                reaction1_b=np.array(reaction1_b)
                reaction1_b_mean=reaction1_b.mean(axis=0)
                reaction1_f=np.array(reaction1_f)
                reaction1_f_mean=reaction1_f.mean(axis=0)



                S_particle=np.array(S_particle)
                #print S_particle.shape
                S_particle_mean=S_particle.mean(axis=0)
                #print S_particle_mean.shape




                ES_particle=np.array(ES_particle)
                ES_particle_mean=ES_particle.mean(axis=0)







                c_s=(S_particle_mean*1.0)/(1.0*self.boxsize**3)
                c_es=(ES_particle_mean*1.0)/(1.0*self.boxsize**3)

                #plt.plot(t,reaction1_b_mean/c_es)
                integral_of_kplus=(reaction1_f_mean*self.boxsize**3)/(c_s*c_es)


                def es_t_modell(t,kplus,kminus,kcomplex):
                    #kminus=how_reaction
                    kcomplex=self.micro_reactionrate_complex
                    kminus=self.micro_reactionrate_backward
                    S_0=(self.particles*1.0)/(self.boxsize**3)
                    E_0=1.0/(self.boxsize**3)
                    #t= np.linspace(0,self.length*self.timestep,10000)
                    S_t=((kminus+kcomplex)/kplus)*lambertw((S_0/((kminus+kcomplex)/kplus))*np.exp((-kminus*E_0*t+S_0)/((kminus+kcomplex)/kplus)))
                    ES_t=((E_0*S_t.real)/(((kminus+kcomplex)/kplus)+S_t.real))*(1-np.exp(-(((kminus+kcomplex)/kplus)+S_t.real)*kplus*t))

                    return ES_t
                best_vals, covar= curve_fit(es_t_modell,t,(ES_particle_mean*1.0)/(1.0*self.boxsize**3))

                ES_model=es_t_modell(t,best_vals[0],1,1)



                #plt.plot(t, (reaction1_f_mean*self.boxsize**3)/(S_particle_mean*ES_particle_mean))

                k1=(np.diff(reaction1_f_mean)*self.boxsize**3)/(S_particle_mean[1:])#*ES_particle_mean[1:])

                #plt.plot(t[1:],k1)

                def makemean(k1,step):
                    i=0
                    index_old=0
                    t_neu=[]
                    t_increment=((step+1)/2)
                    k1_neu=[]
                    while i<self.length/step:
                        index=index_old+step
                        k1_neu.append(k1[index_old:index].mean())
                        t_neu.append(t_increment)
                        t_increment=t_increment+step
                        index_old=index
                        i=i+1
                    return np.array(t_neu), np.array(k1_neu)

                #t_index,k1neu=k1_mean(k1,21)


                #plt.plot(t[t_index],k1neu)
                t_index,k1neu=makemean(k1,101)
                plt.plot(t[t_index],k1neu*20,label=self.alpha)


                #plt.plot(t, (self.boxsize**3)/(ES_particle_mean))
                #plt.plot(t, 1.0/(ES_model))

            elif ii=="kplus":
                reaction1_f=[]

                S_particle=[]
                ES_particle=[]

                for file in os.listdir(self.path+"reaction1/"):
                    f = h5.File(self.path+"reaction1/"+file, "r")
                    reaction1_f.append(f['forwardCounter'][:])
                    t = f['times'][:]
                    f.close()
                for file in os.listdir(self.path+"number_S/"):
                    fn = h5.File(self.path+"number_S/"+file, "r")
                    S_particle.append(fn['particleNumbers'][:])
                    fn.close()
                for file in os.listdir(self.path+"number_C/"):
                    fn = h5.File(self.path+"number_C/"+file, "r")
                    ES_particle.append(fn['particleNumbers'][:])
                    fn.close()


                reaction1_f=np.array(reaction1_f)
                reaction1_f_mean=reaction1_f.mean(axis=0)

                S_particle=np.array(S_particle)
                S_particle_mean=S_particle.mean(axis=0)


                ES_particle=np.array(ES_particle)
                ES_particle_mean=ES_particle.mean(axis=0)

                c_s=(S_particle_mean*1.0)/(1.0*self.boxsize**3)
                c_es=(ES_particle_mean*1.0)/(1.0*self.boxsize**3)

                def s_t_modell(t,kplus,kminus,kcomplex):
                    S_0=(self.particles*1.0)/(self.boxsize**3)
                    E_0=1.0/(self.boxsize**3)
                    kcomplex=self.micro_reactionrate_complex
                    kminus=self.micro_reactionrate_backward
                    S_t=((kminus+kcomplex)/kplus)*lambertw((S_0/((kminus+kcomplex)/kplus))*np.exp((-kcomplex*E_0*t+S_0)/((kminus+kcomplex)/kplus)))
                    return S_t.real
                best_vals, covar= curve_fit(s_t_modell,t,c_s)

                self.k_plus=best_vals[0]

                def es_t_modell(t,kplus,kminus,kcomplex):

                    kcomplex=self.micro_reactionrate_complex
                    kminus=self.micro_reactionrate_backward
                    S_0=(self.particles*1.0)/(self.boxsize**3)
                    E_0=1.0/(self.boxsize**3)
                    S_t=((kminus+kcomplex)/kplus)*lambertw((S_0/((kminus+kcomplex)/kplus))*np.exp((-kminus*E_0*t+S_0)/((kminus+kcomplex)/kplus)))
                    ES_t=((E_0*S_t.real)/(((kminus+kcomplex)/kplus)+S_t.real))*(1-np.exp(-(((kminus+kcomplex)/kplus)+S_t.real)*kplus*t))

                    return ES_t
                best_vals, covar= curve_fit(es_t_modell,t,(ES_particle_mean*1.0)/(1.0*self.boxsize**3))

                ES_model=es_t_modell(t,best_vals[0],1,1)



                #k1=(np.diff(reaction1_f_mean)*self.boxsize**3)/(S_particle_mean[1:])*(1-ES_model[1:]*self.boxsize**3)
                k1=(np.diff(reaction1_f_mean)*self.boxsize**3)/((S_particle_mean[1:])*(1-ES_particle_mean[1:]))

                def makemean(k1,step):
                    i=0
                    index_old=0
                    t_neu=[]
                    t_increment=((step+1)/2)
                    k1_neu=[]
                    while i<self.length/step:
                        index=index_old+step
                        k1_neu.append(k1[index_old:index].mean())
                        t_neu.append(t_increment)
                        t_increment=t_increment+step
                        index_old=index
                        i=i+1
                    return np.array(t_neu), np.array(k1_neu)

                t_index,k1neu=makemean(k1,101)
                tneu=t[t_index]

                def k1_fit(t,h,k_0):
                    if self.alpha==1.0:
                        h=0
                    return k_0*t**(-h)


                h_best, covar_h= curve_fit(k1_fit,tneu[10:],k1neu[10:]*20)
                tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
                for i in range(len(tableau20)):
                    r, g, b = tableau20[i]
                    tableau20[i] = (r / 255., g / 255., b / 255.)
                if self.alpha==1.0:
                    colorthis=tableau20[0]
                if self.alpha==0.9:
                    colorthis=tableau20[1]
                if self.alpha==0.8:
                    colorthis=tableau20[2]
                if self.alpha==0.7:
                    colorthis=tableau20[3]
                if self.alpha==0.6:
                    colorthis=tableau20[4]
                if self.alpha==0.5:
                    colorthis=tableau20[5]
                if self.alpha==0.4:
                    colorthis=tableau20[6]
                if self.alpha==0.3:
                    colorthis=tableau20[7]
                if self.alpha==0.2:
                    colorthis=tableau20[8]
                if self.alpha==0.1:
                    colorthis=tableau20[9]
                if self.alpha==0.05:
                    colorthis=tableau20[10]
                if self.alpha==1.1:
                    colorthis=tableau20[11]
                if self.alpha==1.2:
                    colorthis=tableau20[12]
                if self.alpha==1.3:
                    colorthis=tableau20[13]
                if self.alpha==1.4:
                    colorthis=tableau20[14]
                if self.alpha==1.5:
                    colorthis=tableau20[15]
                if self.alpha==1.6:
                    colorthis=tableau20[16]
                if self.alpha==1.7:
                    colorthis=tableau20[17]
                if self.alpha==1.8:
                    colorthis=tableau20[18]
                if self.alpha==1.9:
                    colorthis=tableau20[19]
                   
                
                
                #tneu=np.array(range(2**16-1))
                #plt.scatter(tneu[np.logspace(0,np.log10(len(tneu)),40).astype('int')[2:]]/self.timestep,p_t_modell(tneu, self.k1_alpha_0,self.micro_reactionrate_backward,self.micro_reactionrate_complex)[np.logspace(0,np.log10(len(tneu)),40).astype('int')[2:]]/S_0,alpha=0.3,color="g",facecolors='none' )

                plt.scatter(t[t_index],k1neu*20,color=colorthis,label="$\\alpha=%.2f$"%self.alpha,facecolors='none')
                plt.loglog(tneu[10:],k1_fit(tneu[10:],h_best[0],h_best[1]),color=colorthis)
                #plt.plot(tneu,np.ones(len(tneu))*self.k_plus,linestyle=":",label="Fit $c_s(t)$")




                print "das ist $h=%.5f$ fuer $\\alpha=%.2f$" %(h_best[0],self.alpha)






            elif ii=="kplus_neu":
                reaction1_f=[]

                S_particle=[]
                ES_particle=[]

                for file in os.listdir(self.path+"reaction1/"):
                    f = h5.File(self.path+"reaction1/"+file, "r")
                    reaction1_f.append(f['forwardCounter'][:])
                    t = f['times'][:]
                    f.close()
                for file in os.listdir(self.path+"number_S/"):
                    fn = h5.File(self.path+"number_S/"+file, "r")
                    S_particle.append(fn['particleNumbers'][:])
                    fn.close()
                for file in os.listdir(self.path+"number_C/"):
                    fn = h5.File(self.path+"number_C/"+file, "r")
                    ES_particle.append(fn['particleNumbers'][:])
                    fn.close()


                reaction1_f=np.array(reaction1_f)
                reaction1_f_mean=reaction1_f.mean(axis=0)

                S_particle=np.array(S_particle)
                S_particle_mean=S_particle.mean(axis=0)


                ES_particle=np.array(ES_particle)
                ES_particle_mean=ES_particle.mean(axis=0)

                c_s=(S_particle_mean*1.0)/(1.0*self.boxsize**3)
                c_es=(ES_particle_mean*1.0)/(1.0*self.boxsize**3)

                def s_t_modell(t,kplus,kminus,kcomplex):
                    S_0=(self.particles*1.0)/(self.boxsize**3)
                    E_0=1.0/(self.boxsize**3)
                    kcomplex=self.micro_reactionrate_complex
                    kminus=self.micro_reactionrate_backward
                    S_t=((kminus+kcomplex)/kplus)*lambertw((S_0/((kminus+kcomplex)/kplus))*np.exp((-kcomplex*E_0*t+S_0)/((kminus+kcomplex)/kplus)))
                    return S_t.real
                best_vals, covar= curve_fit(s_t_modell,t,c_s)

                self.k_plus=best_vals[0]

                def es_t_modell(t,kplus,kminus,kcomplex):

                    kcomplex=self.micro_reactionrate_complex
                    kminus=self.micro_reactionrate_backward
                    S_0=(self.particles*1.0)/(self.boxsize**3)
                    E_0=1.0/(self.boxsize**3)
                    S_t=((kminus+kcomplex)/kplus)*lambertw((S_0/((kminus+kcomplex)/kplus))*np.exp((-kminus*E_0*t+S_0)/((kminus+kcomplex)/kplus)))
                    ES_t=((E_0*S_t.real)/(((kminus+kcomplex)/kplus)+S_t.real))*(1-np.exp(-(((kminus+kcomplex)/kplus)+S_t.real)*kplus*t))

                    return ES_t
                best_vals, covar= curve_fit(es_t_modell,t,(ES_particle_mean*1.0)/(1.0*self.boxsize**3))

                ES_model=es_t_modell(t,best_vals[0],1,1)



                #k1=(np.diff(reaction1_f_mean)*self.boxsize**3)/(S_particle_mean[1:])*(1-ES_model[1:]*self.boxsize**3)
                k1=(np.diff(reaction1_f_mean)*self.boxsize**3)/((S_particle_mean[1:])*(1-ES_particle_mean[1:]))

                def makemean(k1,step):
                    i=0
                    index_old=0
                    t_neu=[]
                    t_increment=((step+1)/2)
                    k1_neu=[]
                    while i<self.length/step:
                        index=index_old+step
                        k1_neu.append(k1[index_old:index].mean())
                        t_neu.append(t_increment)
                        t_increment=t_increment+step
                        index_old=index
                        i=i+1
                    return np.array(t_neu), np.array(k1_neu)

                t_index,k1neu=makemean(k1,301)
                tneu=t[t_index]

                def k1_fit(t,h,k_0):
                    if self.alpha==1.0:
                        h=0
                    return k_0*t**(-h)


                h_best, covar_h= curve_fit(k1_fit,tneu[10:],k1neu[10:]*20)

                k_erban=4.0*np.pi*self.Diffusion*(self.reactiondistance-(np.sqrt(self.Diffusion/(self.micro_reactionrate_forward))*np.tanh(self.reactiondistance*np.sqrt((self.micro_reactionrate_forward)/self.Diffusion))))

                print "Mean kpluscount=",k1neu.mean()*20
                print "std of Mean kpluscount=",k1neu.std()*20/(np.sqrt(len(k1neu)))
                if self.micro_reactionrate_complex==3.0:
                    plt.scatter(t[t_index]/self.timestep,k1neu*20,facecolors='none',color="b",label=" reactions per time interval")
                    #plt.plot(tneu[:]/self.timestep,k1_fit(tneu[:],h_best[0],h_best[1]),linestyle="--",label="mean of reaction counted")
                    plt.semilogy(tneu/self.timestep,np.ones(len(tneu))*self.k_plus,color="b",label="QSSA fitt of $c_s(t)$")
                    plt.plot(tneu/self.timestep,np.ones(len(tneu))*k_erban,linestyle=":",color="b",label="Erban-Chapman")

                else:

                    plt.scatter(t[t_index]/self.timestep,k1neu*20,facecolors='none',color="g")
                    #plt.plot(tneu[:]/self.timestep,k1_fit(tneu[:],h_best[0],h_best[1]),linestyle="--",label="mean of reaction counted")
                    plt.plot(tneu/self.timestep,np.ones(len(tneu))*self.k_plus ,color="g")
                    plt.plot(tneu/self.timestep,np.ones(len(tneu))*k_erban,linestyle=":", color="g")






                print h_best[0],"das ist h"
                print h_best[1],"k_0"
                print self.alpha, "das ist alpha"


            elif ii=="kminus":

                reaction1_b=[]


                ES_particle=[]

                for file in os.listdir(self.path):
                    if  file.endswith("reac_num_S_and_E_to_C.h5"):
                        numberfile1=file[:-24]+"number_c.h5"
                        fc = h5.File(self.path+numberfile1, "r")
                        particles_c = fc['particleNumbers'][:]
                        ES_particle.append(particles_c)
                        f = h5.File(self.path+file, "r")
                        reaction_b=f['backwardCounter'][:]
                        reaction1_b.append(np.float64(reaction_b))
                        t = f['times'][:]
                        f.close()
                        fc.close()

                reaction1_b=np.array(reaction1_b)
                print reaction1_b.shape
                reaction1_b_mean=reaction1_b.mean(axis=0)
                ES_particle=np.array(ES_particle)
                ES_particle_mean=ES_particle.mean(axis=0)
                c_es=(ES_particle_mean*1.0)/(1.0*self.boxsize**3)

                def es_t_modell(t,kplus,kminus,kcomplex):
                    #kminus=how_reaction
                    kcomplex=self.micro_reactionrate_complex
                    kminus=self.micro_reactionrate_backward
                    S_0=(self.particles*1.0)/(self.boxsize**3)
                    E_0=1.0/(self.boxsize**3)
                    #t= np.linspace(0,self.length*self.timestep,10000)
                    S_t=((kminus+kcomplex)/kplus)*lambertw((S_0/((kminus+kcomplex)/kplus))*np.exp((-kminus*E_0*t+S_0)/((kminus+kcomplex)/kplus)))
                    ES_t=((E_0*S_t.real)/(((kminus+kcomplex)/kplus)+S_t.real))*(1-np.exp(-(((kminus+kcomplex)/kplus)+S_t.real)*kplus*t))

                    return ES_t

                best_vals, covar= curve_fit(es_t_modell,t,(ES_particle_mean*1.0)/(1.0*self.boxsize**3))

                ES_model=es_t_modell(t,best_vals[0],1,1)




                kminus=(np.diff(reaction1_b_mean))/(ES_particle_mean[1:])#*ES_particle_mean[1:])


                def makemean(timeseries,step):
                    i=0
                    index_old=0
                    t_neu=[]
                    t_increment=((step+1)/2)
                    timeseries_neu=[]
                    while i<self.length/step:
                        index=index_old+step
                        timeseries_neu.append(timeseries[index_old:index].mean())
                        t_neu.append(t_increment)
                        t_increment=t_increment+step
                        index_old=index
                        i=i+1
                    return np.array(t_neu), np.array(timeseries_neu)

                t_index,kminusneu=makemean(kminus,101)
                plt.plot(t[t_index],kminusneu*20,label=self.alpha)

            elif ii=="kcomplex":

                reaction2_b=[]


                ES_particle=[]

                for file in os.listdir(self.path):
                    if  file.endswith('reac_num_S_and_P_to_C.h5'):
                        numberfile1=file[:-24]+"number_c.h5"
                        fc = h5.File(self.path+numberfile1, "r")
                        particles_c = fc['particleNumbers'][:]
                        ES_particle.append(particles_c)
                        f = h5.File(self.path+file, "r")
                        reaction_b=f['backwardCounter'][:]
                        reaction2_b.append(np.float64(reaction_b))
                        t = f['times'][:]
                        f.close()
                        fc.close()

                reaction2_b=np.array(reaction2_b)
                print reaction2_b.shape
                reaction2_b_mean=reaction2_b.mean(axis=0)
                ES_particle=np.array(ES_particle)
                ES_particle_mean=ES_particle.mean(axis=0)
                c_es=(ES_particle_mean*1.0)/(1.0*self.boxsize**3)

                def es_t_modell(t,kplus,kminus,kcomplex):
                    kcomplex=self.micro_reactionrate_complex
                    kminus=self.micro_reactionrate_backward
                    S_0=(self.particles*1.0)/(self.boxsize**3)
                    E_0=1.0/(self.boxsize**3)
                    S_t=((kminus+kcomplex)/kplus)*lambertw((S_0/((kminus+kcomplex)/kplus))*np.exp((-kminus*E_0*t+S_0)/((kminus+kcomplex)/kplus)))
                    ES_t=((E_0*S_t.real)/(((kminus+kcomplex)/kplus)+S_t.real))*(1-np.exp(-(((kminus+kcomplex)/kplus)+S_t.real)*kplus*t))

                    return ES_t

                best_vals, covar= curve_fit(es_t_modell,t,c_es)
                ES_model=es_t_modell(t,best_vals[0],1,1)




                kcomplex=(np.diff(reaction2_b_mean))/(ES_particle_mean[1:])#*ES_particle_mean[1:])


                def makemean(timeseries,step):
                    i=0
                    index_old=0
                    t_neu=[]
                    t_increment=((step+1)/2)
                    timeseries_neu=[]
                    while i<self.length/step:
                        index=index_old+step
                        timeseries_neu.append(timeseries[index_old:index].mean())
                        t_neu.append(t_increment)
                        t_increment=t_increment+step
                        index_old=index
                        i=i+1
                    return np.array(t_neu), np.array(timeseries_neu)

                t_index,kcomplexneu=makemean(kcomplex,101)
                plt.plot(t[t_index],kcomplexneu*20,label=self.alpha)


            elif ii=="Lineweaver-Burk-Diagramm":
                reaction2_b=[]
                reaction1_f=[]
                particle_number_s=[]
                particle_number_c=[]
                for file in os.listdir(self.path):


                    if  file.endswith('reac_num_S_and_P_to_C.h5'):
                        f = h5.File(self.path+file, "r")
                        reaction_b=f['backwardCounter'][:]
                        t = f['times'][:]
                        reaction2_b.append(np.float64(reaction_b))
                        f.close()
                    elif  file.endswith('reac_num_S_and_E_to_C.h5'):
                        f = h5.File(self.path+file, "r")
                        reaction_f = f['forwardCounter'][:]
                        reaction1_f.append(np.float64(reaction_f))
                        t = f['times'][:]
                        f.close()

                    elif  file.endswith("number_c.h5"):
                        f = h5.File(self.path+file, "r")
                        t = f['times'][:]
                        fc = f['particleNumbers'][:]
                        particle_number_c.append(fc)
                        f.close()


                    elif file.endswith("number_s.h5"):
                        f = h5.File(self.path+file, "r")
                        t = f['times'][:]
                        fc = f['particleNumbers'][:]
                        particle_number_s.append(fc)
                        f.close()

                particle_number_c=np.array(particle_number_c)
                particle_number_c_mean=particle_number_c.mean(axis=0)
                reaction2_b=np.array(np.diff(reaction2_b))
                reaction1_f=np.array(np.diff(reaction1_f))
                reaction1_f_mean=reaction1_f.mean(axis=0)
                particle_number_s=np.array(particle_number_s)
                particle_number_s_mean=particle_number_s.mean(axis=0)
                linestyle='-'
                if self.alpha==0.5:
                    linestyle='--'
                regr = linear_model.LinearRegression()
                X=(self.boxsize**3)/particle_number_s_mean[10:]
                Y=1./(self.micro_reactionrate_complex*particle_number_c_mean[10:])
                X=np.reshape(X,(len(X),1))

                print reaction1_f_mean.shape
                print particle_number_s_mean.shape



    def show_particle(self):
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

            self.observables=["number_of_particle"]
            self.show_observable()
            boxsize1=self.boxsize*1.0
            S_0=(self.particles*1.0)/(boxsize1**3)
            E_0=1.0/(boxsize1**3)
            k2=self.k_minus
            k_1=self.k_plus
            k_complex=self.k_complex
            K_m=(k2+k_complex)/k_1

            """
            S. Schnell and C. Mendoza: Closed form solution for time-dependent enzyme
            kinetics, Journal of Theoretical Biology 187, 207 (1997)
            """
            S_t=K_m*lambertw((S_0/K_m)*np.exp((-k_complex*E_0*self.time+S_0)/K_m))
            ES_t=((E_0*S_t)/(K_m+S_t))*(1-np.exp(-(K_m+S_t)*k_1*self.time))
            E_t=E_0-ES_t
            P_t=S_0-S_t-ES_t
            """
            Plotting
            """
            plt.plot(self.time,ES_t/E_0,linestyle=":" ) #complex
            plt.plot(self.time,E_t/E_0,linestyle=":")  #enzyme
            plt.plot(self.time,S_t/S_0,linestyle=":")  #substrate
            plt.plot(self.time,P_t/S_0,linestyle=":")  #product

    def show_k1(self):
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

            self.observables=["number_of_particle"]
            self.show_observable(show_plot="off")
            boxsize1=self.boxsize*1.0
            S_0=(self.particles*1.0)/(boxsize1**3)
            E_0=1.0/(boxsize1**3)
            k2=self.k_minus
            k_1=self.k_plus #*(complexb.time*0.05)**(np.log(complexb.alpha)/np.pi)
            k_complex=self.k_complex
            K_m=(k2+k_complex)/k_1
            k_erban=4.0*np.pi*self.Diffusion*(self.reactiondistance-(np.sqrt(self.Diffusion/(self.micro_reactionrate_forward))*np.tanh(self.reactiondistance*np.sqrt((self.micro_reactionrate_forward)/self.Diffusion))))
            print "erban chapmann:k1=",k_erban
            print "fitted k1",k_1

    def show_k2(self):
            self.observables=["kminus"]
            self.show_observable()
    def show_particle_alpha(self):
            self.observables=["Michaelis_fit_alpha"]
            self.show_observable()
            return self.k1_alpha_0

    def show_k1_better(self):
            self.observables=["kplus_neu"]
            self.show_observable()
            k_erban=4.0*np.pi*self.Diffusion*(self.reactiondistance-(np.sqrt(self.Diffusion/(self.micro_reactionrate_forward))*np.tanh(self.reactiondistance*np.sqrt((self.micro_reactionrate_forward)/self.Diffusion))))
            print "erban chapmann:k1=",k_erban
            k_1=self.k_plus
            print "fitted k1",k_1


    def show_kcomplex(self):
            self.observables=["kcomplex"]
            self.show_observable()

    def show_radial(self):
            self.observables=["radial"]
            self.show_observable()
    def show_reaction(self):
            self.observables=["reaction"]
            self.show_observable()
    def show_MSD(self):
            self.observables=["MSD"]
            self.show_observable()
            plt.plot(self.time[1:],6*((self.time[1:])**self.alpha)*self.Diffusion,linestyle="--")
