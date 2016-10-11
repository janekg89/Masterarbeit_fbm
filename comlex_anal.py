import logging
import os
import random
import string
import time
#from __future__ import print_function
import h5py as h5
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import complex_sim
from sklearn import datasets, linear_model
from scipy.optimize import curve_fit
from scipy.special import lambertw



import revreaddy as rdy

reload(logging)
logging.basicConfig(
    format='%(asctime)s %(levelname)s: %(message)s',
    datefmt='%Y/%m/%d %I:%M:%S',
    level=logging.DEBUG
    )
class Analy_Complex(complex_sim.Sim_Complex):
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
                        # The function y(x)=a+b*x+c*x^2 is a polynomial
                        # in this example.
                        y=K_m*lambertw((S_0/K_m)*np.exp((-k2*E_0*x+S_0)/K_m))
                        #S_t=S_t.astype(float)
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
                    #k2=self.micro_reactionrate_complex
                    for n in range(len(X)):
                        x = X[n]
                        # The function y(x)=a+b*x+c*x^2 is a polynomial
                        # in this example.
                        y_old=K_m*lambertw((S_0/K_m)*np.exp((-k2*E_0*x+S_0)/K_m))
                        y=(((E_0*y_old)/(K_m+y_old))*(1-np.exp(-(K_m+y_old)*((self.micro_reactionrate_backward+self.micro_reactionrate_complex)/K_m)*x)))
                        #S_t=S_t.astype(float)
                        chi2 = chi2 + (Y[n] - y)**2
                    return chi2
    def func_p(self,params, X, Y):
                    # extract current values of fit parameters from input array
                    K_m = params[0]
                    #k_1=params[1]
                    # compute chi-square
                    chi2 = 0.0
                    S_0=(self.particles*1.0)/(1.0*self.boxsize**3)
                    E_0=1.0/(1.0*self.boxsize**3)
                    k2=self.micro_reactionrate_complex
                    for n in range(len(X)):
                        x = X[n]
                        # The function y(x)=a+b*x+c*x^2 is a polynomial
                        # in this example.
                        y_old=K_m*lambertw((S_0/K_m)*np.exp((-k2*E_0*x+S_0)/K_m))
                        y=(((E_0*y_old)/(K_m+y_old))*(1-np.exp(-(K_m+y_old)*((self.micro_reactionrate_backward+self.micro_reactionrate_complex)/K_m)*x)))
                        y_neu=S_0-y_old-y
                        #S_t=S_t.astype(float)
                        chi2 = chi2 + (Y[n] - y_neu)**2


                    return chi2
    def show_observable(self):
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
                color="y"
                if self.alpha==0.5:
                    color="b"
                    linestyle='--'
                    marker=""
                    #markersize=1.0

                #plt.plot(t, particle_number_c_mean,marker=marker,linestyle=linestyle,color=color,markevery=every ,label="complex") #of complex

                particle_number_s=np.array(particle_number_s)
                particle_number_s_mean=particle_number_s.mean(axis=0)
                plt.plot(t, particle_number_s_mean/(self.particles*1.0),linestyle=linestyle,color="r" ,marker=marker,markevery=every, label="substrate") #of substrate

                particle_number_p=np.array(particle_number_p)
                particle_number_p_mean=particle_number_p.mean(axis=0)

                #v=np.diff(particle_number_p*1.0)
                #v_mean=v.mean(axis=0)
                #print self.micro_reactionrate_complex
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
                #print covar1[1]





                self.k_all=[best_vals[0],best_vals1[0],best_vals2[0]]
                self.var_all=[covar[0],covar1[0],covar2[0]]

                #self.k_plus=(best_vals2[0]+best_vals[0]+best_vals1[0])/3.0
                self.k_plus=best_vals1[0]
                self.k_minus=self.micro_reactionrate_backward
                self.k_complex=self.micro_reactionrate_complex



                from scipy.optimize import fmin as simplex
                # "fmin" is not a sensible name for an optimisation package.
                # Rename fmin to "simplex"

                # Define the objective function to be minimised by Simplex.
                # params ... array holding the values of the fit parameters.
                # X      ... array holding x-positions of observed data.
                # Y      ... array holding y-values of observed data.
                # Err    ... array holding errors of observed data.





                xdata = t*1.0
                ydata = (particle_number_s_mean*1.0)/self.boxsize**3
                ydata_es = (particle_number_c_mean*1.0)/self.boxsize**3
                ydata_p = (particle_number_p_mean*1.0)/self.boxsize**3


                #Initial guess.
                x0    = [1.0]

                # Apply downhill Simplex algorithm.
                #self.K_m_mean = simplex(self.func, x0, args=(xdata, ydata))
                #print  simplex(self.func_es, x0, args=(xdata, ydata_es))
                #self.K_m_mean = simplex(self.func_p, x0, args=(xdata, ydata_p))

                #print self.K_m_mean_1-self.K_m_mean

                #self.k_plus=(self.k2+self.k2)/self.K_m_mean

                #print self.func(self.K_m_mean,t,particle_number_s_mean/self.boxsize**3)
                #print self.func(0.5*self.K_m_mean,t,particle_number_s_mean/self.boxsize**3)








                #plt.plot(1.0*particle_number_s_mean[1:]/(self.boxsize**3),v_mean,label="Michaelis Menten Diagramm",marker=".", linestyle="")



                #b = 0.5-0.5/np.exp(self.micro_reactionrate_forward*(1-particle_number_c_mean.max())*(t))
                #substrate=1*np.exp(-self.micro_reactionrate_backward*(1-particle_number_c_mean)/0.59)
                #d=(particle_number_s_mean/(self.particles*0.59))/(1+particle_number_s_mean/(self.particles *0.59))
                #plt.plot(t,substrate,label="coulomb")


                #plt.plot(t, particle_number_p_mean/(self.particles*1.0),linestyle=linestyle,color="g"  ,marker=marker,markevery=every,label=" product " ) # of product

                #plt.plot(t, particle_number_s_mean/boxsize**3,linestyle=linestyle,label=" michaelis " )

                #k_m=((particle_number_s_mean)/boxsize**3)*((1.0/(particle_number_c_mean*1.0)-1))

                #self.K_m_mean= k_m[100:200].mean()

                #self.k_plus=(self.micro_reactionrate_backward+self.micro_reactionrate_complex)/self.K_m_mean


                #plt.plot(t, k_m,linestyle=linestyle,label=" michaelis " )




                #plt.plot(t, np.cumsum(particle_number_c_mean*self.micro_reactionrate_complex)/(self.particles*20.0),linestyle=":",label=" boxsize = %.2f " %self.boxsize)

                #plt.plot(t[:-1], np.diff(particle_number_p_mean)/self.timestep,label="product")

            elif ii== "radial":
                radial_s=[]
                radial_p=[]
                particle_number_s=[]
                particle_number_p=[]
                linestyle='-'
                if self.alpha==0.5:
                    linestyle='--'
                for file in os.listdir(self.path):
                    if  file.endswith("radialS.h5"):
                        numberfile=file[:-10]+"number_s.h5"
                        fn = h5.File(self.path+numberfile, "r")

                        numberfile1=file[:-10]+"number_p.h5"
                        fn1 = h5.File(self.path+numberfile1, "r")
                        particle_number_p.append(fn1["particleNumbers"][-1])
                        #print numberfile

                        f = h5.File(self.path+file, "r")
                        bins = f['bins'][:]
                        bincenters = f['binCenters'][:]

                        if not bins.sum()==0:
                            #print fn["particleNumbers"][-1]
                            radial_s.append(bins)
                            #radial_s.append(bins/bins.sum())
                        particle_number_s.append(fn["particleNumbers"][-1])
                        #else:

                            #print bins.sum()

                        f.close()
                        fn.close()
                        fn1.close()
                    elif  file.endswith("radialP.h5"):
                        f = h5.File(self.path+file, "r")
                        #print f.keys()
                        bins = f['bins'][:]
                        #print bins.sum()
                        bincenters = f['binCenters'][:]
                        #print bincenters
                        if not bins.sum()==0:
                            #print fn["particleNumbers"][-1]
                            #radial_p.append(bins/bins.sum())
                            radial_p.append(bins)
                        #else:
                            #print bins.sum()
                        f.close()



                particle_number_s=np.array(particle_number_s)
                particle_number_p=np.array(particle_number_p)
                #if (particle_number_s+particle_number_p).any() is not 20:
                #    print particle_number_s+particle_number_p


                #print particle_number_s+particle_number_p
                radial_s=np.array(radial_s)
                #print radial_s.mean(axis=0).shape
                plt.plot(bincenters,(radial_s.mean(axis=0)*30.0)/radial_s.mean(axis=0).sum())
                radial_p=np.array(radial_p)

                plt.plot(bincenters,(radial_p.mean(axis=0)*30.0)/radial_p.mean(axis=0).sum(),linestyle="--")
                #radial_s_norm=(radial_s*1.0)/(particle_number_s[:, None])
                #print particle_number_s[:, None].min()
                #print radial_s_norm.max()
                #print radial_s.max(axis=1)
                #print np.where(radial_s == radial_s.max())

                #print radial_s.sum(axis=1)


                #radial_s_norm=(radial_s*1.0)/(self.particles*1.0)
                #print radial_s_norm.shape
                #radial_s_norm_mean=radial_s_norm.mean(axis=0)
                #print radial_s_norm_mean.shape






                #radial_s_mean=radial_s.mean(axis=0)





                #radialfunction_s = ((self.boxsize**3)*radial_s_norm_mean)/((self.boxsize/2.0-self.reactiondistance)/(50.0)*8*np.pi)# 8 because of checking for radial distribution next to complex and substrate
                #if len(radial_s) > 0:
                    #plt.plot(bincenters,radialfunction_s ,linestyle=linestyle,label="substrate")


                #radial_p_norm=radial_p/(particle_number_p[:, None])
                #radial_p_norm_mean=radial_p_norm.mean(axis=0)



                #radial_p_mean=radial_p.mean(axis=0)

                #radialfunction_p = (self.boxsize**3*radial_p_norm_mean)/((self.boxsize/2.0-self.reactiondistance)/(50.0)*8*np.pi) # 8 because of checking for radial distribution next to complex and substrate
                #if len(radial_p) > 0:
                    #plt.plot(bincenters, radialfunction_p,linestyle=linestyle,label="product")

            elif ii=="MSD":
                msd_s=[]

                for file in os.listdir(self.path):
                    if  file.endswith("msd.h5"):
                        f = h5.File(self.path+file, "r")
                        msd = f['meanSquaredDisplacements'][:]
                        t = f['times'][:]
                        if len(msd)-1 == self.length:
                            msd_s.append(np.array(msd))

                        f.close()
                msd_s=np.array(msd_s)
                msd_s_mean=msd_s.mean(axis=0)
                plt.plot(t[1:],msd_s_mean[1:],label="substrate")

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
                #print t[-1]
               # print (reaction1_b_mean/c_es)[-1]/t[-1]
                #print c_es


                #plt.plot(t,reaction1_f_mean,linestyle="--",label=self.micro_reactionrate_backward)
                #plt.plot(t,reaction2_b_mean,label=self.micro_reactionrate_backward)

                #plt.plot(t,S_particle_mean,label=self.micro_reactionrate_backward)
                #plt.plot(t,reaction1_b_mean+reaction2_b_mean)



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
                best_vals, covar= curve_fit(es_t_modell,t,(ES_particle_mean*1.0)/(1.0*self.boxsize**3))

                ES_model=es_t_modell(t,best_vals[0],1,1)



                #plt.plot(t, (reaction1_f_mean*self.boxsize**3)/(S_particle_mean*ES_particle_mean))

                k1=np.diff(reaction1_f_mean)/(S_particle_mean[1:]*ES_model[1:])

                plt.plot(t[1:],k1)

                def k1_mean(k1,step):
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

                t_index,k1neu=k1_mean(k1,21)


                plt.plot(t[t_index],k1neu)
                t_index,k1neu=k1_mean(k1,101)
                plt.plot(t[t_index],k1neu)


                #plt.plot(t, (self.boxsize**3)/(ES_particle_mean))
                #plt.plot(t, 1.0/(ES_model))





















                """
                k_mean=[]
                print ES_particle.shape
                for i in range(len(ES_particle[0,:])):
                    iterator=0
                    k_sum=0
                    for ii in range(len(ES_particle[:,0])):

                        if ES_particle[ii,i]==0:
                            k_sum=k_sum+((reaction1_f[ii,i]*self.boxsize**3)/(S_particle[ii,i]))
                            iterator=iterator+1
                    k_mean.append(k_sum/iterator)

                k_mean=np.array(k_mean)

                def kplus_modell(t,kplus,a):
                    a=0
                    return (a+t*kplus)*ES_particle.mean()


                best_vals, covar= curve_fit(kplus_modell,t,k_mean*1.0)

                print best_vals,covar



                #plt.plot(t,k_mean*ES_particle.mean())
                #plt.plot(t,ES_particle.mean()*kplus_modell(t,best_vals[0],best_vals[1]))


                #ES_particle_mean=ES_particle.mean(axis=0)


                #plt.plot(np.diff((reaction1_f_mean*self.boxsize**3)/(S_particle_mean*ES_particle_mean))*20)




                assert (len(reaction1_f_mean)!=0),"eather no simulation or simulation observable was not ;all;!"


                def kplus_modell(t,kplus,a):
                    a=0
                    return a+t*kplus


                #best_vals, covar= curve_fit(kplus_modell,t,reaction1_f_mean*1.0)

                #print best_vals,covar




                '''
                for i in range(int(self.length*self.timestep)):


                    index=i*20
                    print index
                    reaction1_f_new.append(reaction1_f[:,index+1]-reaction1_f[:,index])


                reaction1_f_new=(np.diff(np.array(reaction1_f_new),axis=1)).mean(axis=1)
                '''

                #plt .plot(t[1:],reaction1_f_mean)
               # plt.plot(t,kplus_modell(t,best_vals[0],best_vals[1]))

                '''
                reaction1_f_mean=reaction1_f.mean(axis=0)
                newreactiondiff=[]
                indexold=0
                print int(self.timestep*self.length)
                for i in range(int(self.timestep*self.length)):
                    index=int(i/self.timestep)
                    #print index
                    newreactiondiff.append(-reaction1_f[:,indexold]+reaction1_f[:,index])
                    indexold=index

                #reaction1diff_f_mean=

                #plt.plot(range(int(self.timestep*self.length)),np.array(newreactiondiff).mean(axis=1))
                #print reaction1_f_mean
                '''
                linestyle='-'
                if self.alpha==0.5:
                    linestyle='--'
                #K=(reaction2_b_mean+reaction1_b_mean)/reaction1_f_mean

                #plt.plot(t,reaction1_f_mean/self.timestep,linestyle=linestyle,label="forward_1 a=%.2f " %self.alpha)
                #plt.plot(t[:],reaction1_b_mean,linestyle=linestyle,label="backward_1 a=%.2f " %self.alpha)
                #plt.plot(t[:],reaction2_b_mean,linestyle=linestyle,label="backward_2 a=%.2f " %self.alpha)
                #plt.plot(t[:],K,linestyle=linestyle,label="Michaelis Konstante")
                #plt.plot(t[:],reaction1_f_mean/reaction1_b_mean,linestyle=linestyle,label="ratio forwad backward a=%.2f " %self.alpha) #ratio forward backward

                """

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
                reaction2_b_mean=reaction2_b.mean(axis=0)

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

                #regr.fit(X,Y)
                #print('V_max', 1./regr.intercept_)

                #plt.plot(t[1:],((particle_number_s_mean[1:]*1.0)/(1.0*self.boxsize**3))*(1.0*reaction1_f_mean/self.timestep), label="konstant")



                #plt.plot((self.boxsize**3)/particle_number_s_mean,1./(self.micro_reactionrate_complex*particle_number_c_mean),linestyle=linestyle,label="a=%.2f " %self.alpha)
                #plt.plot(X,regr.predict(X),linestyle=linestyle,label="a=%.2f " %self.alpha)













particlenumber=150
boxsize=7.0
D=1.0/6.0
R=1.0
alpha=0.5
lambda_plus=1.0
lambda_c=1.0
lambda_minus=1.0
length=2**13





'''
for particlenumber in [10,50,100,500,1000]:
    boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/200)+(4.0/3.0*np.pi))**(1/3.0)
    complex=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.5,length,boxsize,particlenumber,0.05)
    complex.observables=["number_of_particle","radial","MSD","reaction"]
    for index1 in range(500):
        print R
        print index1
        complex.run()

    complex1=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,1.0,length,boxsize,particlenumber,0.05)
    complex1.observables=["number_of_particle","radial","reaction","MSD"]
    for index2 in range(500):
        print index2
        complex1.run()


n=0
for particlenumber in [50,100,500,1000]:
    n=+1
    print n
    boxsize=(((7**3-((4./3.0)*np.pi))*particlenumber/200)+(4.0/3.0*np.pi))**(1/3.0)
    complexa=Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.5,length,boxsize,particlenumber,0.05)
    complexa.observables=["Lineweaver-Burk-Diagramm"]
    complexa.show_observable()

    complexb=Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,1.0,length,boxsize,particlenumber,0.05)
    complexb.observables=["Lineweaver-Burk-Diagramm"]
    complexb.show_observable()

'''

'''
complex=complex_sim.Sim_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.7,2**8,boxsize,particlenumber,0.05)
complex.observables=["radial"]
for index1 in range(10):
        print index1
        complex.run()

complexb=Analy_Complex(D,R,lambda_plus,lambda_minus,lambda_c,0.7,2**8,boxsize,particlenumber,0.05)
complexb.observables=["radial"]
complexb.show_observable()

plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()
'''

