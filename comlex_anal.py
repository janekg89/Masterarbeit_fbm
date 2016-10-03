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
                if self.alpha==0.5:
                    linestyle='--'

                plt.plot(t, particle_number_c_mean,linestyle=linestyle,label="complex") #of complex

                particle_number_s=np.array(particle_number_s)
                particle_number_s_mean=particle_number_s.mean(axis=0)
                plt.plot(t, particle_number_s_mean/(self.particles*1.0),label="substrate") #of substrate

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


                plt.plot(t, particle_number_p_mean/(self.particles*1.0),linestyle=linestyle,label=" product " ) # of product

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
                        #print numberfile

                        f = h5.File(self.path+file, "r")
                        bins = f['bins'][:]
                        bincenters = f['binCenters'][:]
                        if not fn["particleNumbers"][-1]==0:
                            radial_s.append(bins)
                            particle_number_s.append(fn["particleNumbers"][-1])

                        f.close()
                        fn.close()
                    elif  file.endswith("radialP.h5"):
                        numberfile=file[:-10]+"number_p.h5"
                        fn = h5.File(self.path+numberfile, "r")


                        f = h5.File(self.path+file, "r")
                        #print f.keys()
                        bins = f['bins'][:]
                        bincenters = f['binCenters'][:]
                        #print bincenters
                        if not fn["particleNumbers"][-1]==0:
                            particle_number_p.append(fn["particleNumbers"][-1])
                            radial_p.append(bins)
                        f.close()
                        fn.close()

                particle_number_s=np.array(particle_number_s)
                particle_number_p=np.array(particle_number_p)





                radial_s=np.array(radial_s)


                radial_s_norm=(radial_s*1.0)/(particle_number_s[:, None])
                #radial_s_norm=(radial_s*1.0)/(self.particles*1.0)
                #print radial_s_norm.shape
                radial_s_norm_mean=radial_s_norm.mean(axis=0)
                #print radial_s_norm_mean.shape






                #radial_s_mean=radial_s.mean(axis=0)





                radialfunction_s = ((self.boxsize**3)*radial_s_norm_mean)/((self.boxsize/2.0-self.reactiondistance)/(50.0)*8*np.pi)# 8 because of checking for radial distribution next to complex and substrate
                if len(radial_s) > 0:
                    plt.plot(bincenters,radialfunction_s ,label="substrate")

                radial_p=np.array(radial_p)
                radial_p_norm=radial_p/(particle_number_p[:, None])
                radial_p_norm_mean=radial_p_norm.mean(axis=0)



                #radial_p_mean=radial_p.mean(axis=0)

                radialfunction_p = (self.boxsize**3*radial_p_norm_mean)/((self.boxsize/2.0-self.reactiondistance)/(50.0)*8*np.pi) # 8 because of checking for radial distribution next to complex and substrate
                if len(radial_p) > 0:
                    plt.plot(bincenters, radialfunction_p,linestyle=linestyle,label="product")

            elif ii=="MSD":
                msd_s=[]
                for file in os.listdir(self.path):
                    if  file.endswith("msd.h5"):
                        f = h5.File(self.path+file, "r")
                        msd = f['meanSquaredDisplacements'][:]
                        t = f['times'][:]
                        msd_s.append(msd)
                        f.close()
                msd_s=np.array(msd_s)
                msd_s_mean=msd_s.mean(axis=0)
                plt.plot(t,msd_s_mean,label="substrate")


            elif ii=="reaction":
                reaction1_f=[]
                reaction1_b=[]
                reaction2_b=[]
                for file in os.listdir(self.path):
                    if  file.endswith("reac_num_S_and_E_to_C.h5"):
                        f = h5.File(self.path+file, "r")
                        reaction_f = f['forwardCounter'][:]
                        reaction_b=f['backwardCounter'][:]
                        reaction1_b.append(np.float64(reaction_b))
                        reaction1_f.append(np.float64(reaction_f))
                        t = f['times'][:]
                        f.close()


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
                #print reaction1_f_mean
                linestyle='-'
                if self.alpha==0.5:
                    linestyle='--'
                K=(reaction2_b_mean+reaction1_b_mean)/reaction1_f_mean
                #plt.plot(t[:],reaction1_f_mean,linestyle=linestyle,label="forward_1 a=%.2f " %self.alpha)
                #plt.plot(t[:],reaction1_b_mean,linestyle=linestyle,label="backward_1 a=%.2f " %self.alpha)
                plt.plot(t[:],reaction2_b_mean,linestyle=linestyle,label="backward_2 a=%.2f " %self.alpha)
                #plt.plot(t[:],K,linestyle=linestyle,label="Michaelis Konstante")
                plt.plot(t[:],reaction1_f_mean/reaction1_b_mean,linestyle=linestyle,label="ratio forwad backward a=%.2f " %self.alpha) #ratio forward backward



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

