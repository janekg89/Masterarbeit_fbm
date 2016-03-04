# -*- coding: utf-8 -*-
__author__ = 'janek'
import numpy as np # module for scientific computing
from scipy import integrate
import matplotlib.pyplot as plt # module for plotting "a la" matlab

class Felix_Method():

    def __init__(self, D, particles,length,alpha,dt):
        self.D=D
        self.K_alpha=D
        self.particles=particles
        self.n=length
        self.alpha=alpha
        self.ki=self.n*2
        self.frq=np.fft.fftfreq(self.ki)*(np.pi*2./(dt))
        self.t=np.array(range(length))
        self.dt=dt

    def z(self):
        """
        :param D: Diffusionskoeffizient

        :param alpha: Anomalieparameter



        :return: z: Correlationsfunktion im Frequenzraum welche auf die gewöhnliche diffusion angewendet wird.


        """
        z=((((+1j*self.frq))**(1-self.alpha))*self.K_alpha*np.math.gamma(1+self.alpha))*np.exp(np.pi/self.ki)
        return z

    def compute_trajectory(self):
        """
        :param D: Diffusionskoeffizient
        :param particles: Anzahl an Teilechen welche simuliert werden sollen.

        :param length: Länge der Trajektorien, welche simuliert werden sollen.

        :param alpha: Anomalieparameter

        :return: ( Particles x Length ) Array of Trajectories of all particles. Close to :cite:`Graigmile2003`


        """
        r_t_allparticles=[]
        #plt.plot(self.frq,"r")
        #plt.plot( np.arange(-np.pi/self.dt , np.pi/self.dt,(np.pi/(self.n*self.dt))))
        #plt.show()
        for particle in range(self.particles):
            v_t=(np.random.normal(0,np.sqrt(self.dt),size=self.ki)) #todo mit matrix.shape können man eine zufällige verteilung von allen daten generieren
            v_t=np.array(v_t)
            v_frq=np.fft.fft(v_t)
            v_ano_frq= np.sqrt(self.z().real*2.)*v_frq
            v_ano_frq[0]=np.random.normal(0,np.sqrt(2.*self.K_alpha*(self.ki*self.dt)**self.alpha))
            #v_ano_frq[self.n]=np.sqrt(self.z()[self.n].real*self.n*2)*v_t[self.ki-1]
            #v_ano_frq[self.n-1]=np.sqrt(self.z()[self.n].real*self.n*2)*v_t[self.ki-1]

            v_ano_t=np.fft.ifft(v_ano_frq)
            #r_t1=np.cumsum(v_ano_t[:self.n].real) #Ort bei anomaler Diffusion in Abhängigkeit von der zeit
            r_t=integrate.cumtrapz(v_ano_t[:self.n].real,initial=0) #Ort bei anomaler Diffusion in Abhängigkeit von der zeit
            r_t_allparticles.append(r_t) # Trajektorie bei anomaler Diffusion für alle teilchen
        return np.array(r_t_allparticles)





