# -*- coding: utf-8 -*-
__author__ = 'janek'
import numpy as np # module for scientific computing
import matplotlib.pyplot as plt # module for plotting "a la" matlab

class Felix_Method():
    """
    Simulation von Trajektorien der Länge n von mit einer Anzahl (particles) an Teilchen von Anomaler Diffusion.
    Die Methode benutzt eine FFT funktioniert ca. für alpha >0.25
    """
    def __init__(self,D,particles,length,alpha):
        self.D=D
        self.K_alpha=D
        self.particles=particles
        self.n=length
        self.alpha=alpha
        self.ki=self.n*2
        self.frq=np.fft.fftfreq(self.ki)*np.pi*2
        self.t=np.array(range(length))

    def z(self):
        """
        die Correlationsfunktion im Frequenzraum welche auf die gewöhnliche diffusion angewendet wird
        """
        z=(((+1j*self.frq)**(1-self.alpha))*self.K_alpha*np.math.gamma(1+self.alpha))
        return z

    def msdanalyt(self):
        return 2*((self.t)**self.alpha)*self.K_alpha

    def compute_trajectory(self):
        r_t_allparticles=[]
        for particle in range(self.particles):
            v_t=(np.random.normal(0,1,size=self.ki))
            v_t=np.array(v_t)
            v_frq=np.fft.fft(v_t)
            v_ano_frq= np.sqrt(Felix_Method(self.D,self.particles,self.n,self.alpha).z().real*2.)*v_frq
            v_ano_frq[0]=np.random.normal(0,np.sqrt(2.*self.K_alpha*self.ki**self.alpha))
            v_ano_frq[self.n]=np.sqrt(Felix_Method(self.D,self.particles,self.n,self.alpha).z()[self.n].real*self.n*2)*v_t[self.ki-1]
            #v_ano_frq[self.n-1]=np.sqrt(Felix_Method(self.D,self.particles,self.n-1,self.alpha).z()[self.n].real*self.n*2)*v_t[self.ki-1]

            v_ano_t=np.fft.ifft(v_ano_frq)
            r_t=np.cumsum(v_ano_t[:self.n].real) #Ort bei anomaler Diffusion in Abhängigkeit von der zeit
            r_t_allparticles.append(r_t) # Trajektorie bei anomaler Diffusion für alle teilchen
        return r_t_allparticles

    def msd_ensemble(self):
        r_t_allparticles=Felix_Method(self.D,self.particles,self.n,self.alpha).compute_trajectory()
        r_t_allparticles=np.array(r_t_allparticles)
        r_t_allparticles_squared=abs(r_t_allparticles)**2
        msd=r_t_allparticles_squared.mean(axis=0)
        std=r_t_allparticles.std(axis=0)
        return  msd, std

    def plotting(self,error=2):
        msd_ensemble,std=Felix_Method(self.D,self.particles,self.n,self.alpha).msd_ensemble()
        colors=['r','b','g','k','c','w','b','r','g','b','k','c','w','b','r','g','b','k','c','w','bo','ro','go','bo','ko','co','wo','bo']
        #fig=plt.plot(range(msd_ensemble.size), msd_ensemble ,colors[2], label="ensemble msd")
        plt.loglog(self.t,Felix_Method(self.D,self.particles,self.n,self.alpha).msdanalyt(),":",color=colors[1], label="analytisch D=%f,particles=%d,length=%d,alpha=%f" %(self.D,self.particles,self.n,self.alpha))
        fig=plt.errorbar(range(msd_ensemble.size), msd_ensemble, yerr=error*std,label="Spektrale Methode mit D=%f,particles=%d, length=%d ,alpha=%f, Std=%f" %(self.D,self.particles,self.n,self.alpha,error))
        plt.legend(loc=2)
        plt.xlabel('Steps', fontsize=14)
        plt.ylabel('MSD', fontsize=14)
        return fig





