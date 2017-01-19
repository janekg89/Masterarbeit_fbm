# -*- coding: utf-8 -*-
__author__ = 'janek'
import numpy as np # module for scientific computing
from scipy import integrate
import matplotlib.pyplot as plt # module for plotting "a la" matlab
import test_cython.genereatefracincrements as ginc
import lown_cython.genereatefracincrements as ginc1
import hatre_cyton.genereatefracincrements as ginc3
import hatre_cyton1.genereatefracincrements as ginc4




class Felix_Method():


    def __init__(self, D, particles,length,alpha,dt,x=2, version="python"):
        self.D=D
        self.K_alpha=D
        self.particles=particles
        self.n=length
        self.alpha=alpha
        self.ki=self.n*x
        self.frq=np.fft.fftfreq(self.ki)*(np.pi*2./(dt))
        self.frq1=np.fft.fftfreq(self.n)*(np.pi*2./(dt))

        self.t=np.array(range(length))
        self.dt=dt
        self.version=version

    def z(self):
        """
        :param D: Diffusionskoeffizient

        :param alpha: Anomalieparameter



        :return: z: Correlationsfunktion im Frequenzraum welche auf die gewöhnliche diffusion angewendet wird.


        """
        z=((((+1j*self.frq))**(1-self.alpha))*self.K_alpha*np.math.gamma(1+self.alpha))
        return z
    def z1(self):
        """
        :param D: Diffusionskoeffizient

        :param alpha: Anomalieparameter



        :return: z: Correlationsfunktion im Frequenzraum welche auf die gewöhnliche diffusion angewendet wird.


        """
        z=((((+1j*self.frq1))**(1-self.alpha))*self.K_alpha*np.math.gamma(1+self.alpha))
        return z

    def z_alternativ(self):
        MSD=(self.t*self.dt)**(self.alpha)*self.K_alpha*2
        # hier kommt die die zweite ableitung der MSD rein
        d2MSD=self.alpha*(self.alpha-1)*(np.array(range(self.ki))*self.dt)**(self.alpha-2.)*self.K_alpha
        z_omega=np.fft.fft(d2MSD[:])
        return z_omega
    def z_alterantive_2(self):
        cov=[]
        for i in self.t[1:-2]:
            cov.append(((i-1)**self.alpha-2*i**self.alpha+(i+1)**self.alpha)/(2))
        cov=np.array(cov)
        #return np.cumsum(np.fft.fft(cov))
        return cov
    def compute_trajectory(self):
        """
        :param D: Diffusionskoeffizient
        :param particles: Anzahl an Teilechen welche simuliert werden sollen.

        :param length: Länge der Trajektorien, welche simuliert werden sollen.

        :param alpha: Anomalieparameter

        :return: ( Particles x Length ) Array of Trajectories of all particles. Close to :cite:`Graigmile2003`


        """

        #plt.plot(self.frq,"r")
        #plt.plot( np.arange(-np.pi/self.dt , np.pi/self.dt,(np.pi/(self.n*self.dt))))
        #plt.show()
        r_t_allparticles=[]
        if self.version == "cpp":
            inc = ginc.pyIncrements(self.n,self.particles)
            inc.generateIncrements(self.D, self.dt, self.alpha)
            a =inc.returnIncrements()
            r_t=np.zeros((self.particles,self.n))
            r_t[:,1:]=np.cumsum(a[:,0,0:self.n-1],axis=1)
            return r_t
        if self.version == "python":
            sqrt2zreal=np.sqrt(self.z().real*2.)
            for particle in range(self.particles):
                #r = np.random.RandomState(1234)
                v_t=(np.random.normal(0,np.sqrt(self.dt),size=self.ki)) #todo mit matrix.shape könnte man eine zufällige verteilung von allen daten generieren
                v_t=np.array(v_t)
                v_frq=np.fft.fft(v_t)
                v_ano_frq= sqrt2zreal*v_frq
                v_ano_frq[0]=np.random.normal(0,np.sqrt(2.*self.K_alpha*(self.ki*self.dt)**self.alpha))
                #v_ano_frq[self.n]=np.sqrt(self.z()[self.n].real*self.n*2)*v_t[self.ki-1]
                #v_ano_frq[self.n-1]=np.sqrt(self.z()[self.n].real*self.n*2)*v_t[self.ki-1]
                v_ano_t=np.fft.ifft(v_ano_frq)
                #r_t1=np.cumsum(v_ano_t[:self.n].real) #Ort bei anomaler Diffusion in Abhängigkeit von der zeit
                r_t=np.zeros(self.n)
                r_t[1:]=np.cumsum(v_ano_t[:self.n].real)[:self.n-1] #Ort bei anomaler Diffusion in Abhängigkeit von der zeit
                r_t_allparticles.append(r_t) # Trajektorie bei anomaler Diffusion für alle teilchen
                return np.array(r_t_allparticles)
        if self.version == "python_notrick":
            sqrt2zreal=np.sqrt(self.z1().real*2.)
            for particle in range(self.particles):
                #r = np.random.RandomState(1234)
                v_t=(np.random.normal(0,np.sqrt(self.dt),size=self.n)) #todo mit matrix.shape könnte man eine zufällige verteilung von allen daten generieren
                v_t=np.array(v_t)
                v_frq=np.fft.fft(v_t)
                v_ano_frq= sqrt2zreal*v_frq
                #v_ano_frq[0]=np.random.normal(0,np.sqrt(2.*self.K_alpha*(self.ki*self.dt)**self.alpha))
                #v_ano_frq[self.n]=np.sqrt(self.z()[self.n].real*self.n*2)*v_t[self.n-1]
                #v_ano_frq[self.n-1]=np.sqrt(self.z()[self.n].real*self.n*2)*v_t[self.n-1]
                v_ano_t=np.fft.ifft(v_ano_frq)
                #r_t1=np.cumsum(v_ano_t[:self.n].real) #Ort bei anomaler Diffusion in Abhängigkeit von der zeit
                r_t=np.zeros(self.n)
                r_t[1:]=np.cumsum(v_ano_t[:self.n].real)[:self.n-1] #Ort bei anomaler Diffusion in Abhängigkeit von der zeit
                r_t_allparticles.append(r_t) # Trajektorie bei anomaler Diffusion für alle teilchen
            return np.array(r_t_allparticles)
        if self.version == "python_trick":
            sqrt2zreal=np.sqrt(self.z1().real*2.)
            for particle in range(self.particles):
                #r = np.random.RandomState(1234)
                v_t=(np.random.normal(0,np.sqrt(self.dt),size=self.n)) #todo mit matrix.shape könnte man eine zufällige verteilung von allen daten generieren
                v_t=np.array(v_t)
                v_frq=np.fft.fft(v_t)
                v_ano_frq= sqrt2zreal*v_frq
                v_ano_frq[0]=np.random.normal(0,np.sqrt(2.*self.K_alpha*(self.ki*self.dt)**self.alpha))
                #v_ano_frq[self.n]=np.sqrt(self.z()[self.n].real*self.n*2)*v_t[self.ki-1]
                #v_ano_frq[self.n-1]=np.sqrt(self.z()[self.n].real*self.n*2)*v_t[self.n-1]
                v_ano_t=np.fft.ifft(v_ano_frq)
                #r_t1=np.cumsum(v_ano_t[:self.n].real) #Ort bei anomaler Diffusion in Abhängigkeit von der zeit
                r_t=np.zeros(self.n)
                r_t[1:]=np.cumsum(v_ano_t[:self.n].real)[:self.n-1] #Ort bei anomaler Diffusion in Abhängigkeit von der zeit
                r_t_allparticles.append(r_t) # Trajektorie bei anomaler Diffusion für alle teilchen
            return np.array(r_t_allparticles)
        if self.version=="lowen":
            r_t_allparticles=[]
            n=np.array(range(2*self.n))*1.
            r_x1=(self.dt**self.alpha)*(1-(n[:self.n+1]/self.n)**self.alpha)/2
            r_x2=r_x1[::-1]

            r_x=np.append(r_x1,r_x2[1:])
            s_x=np.fft.fft(r_x)
            X_k1=0.0+1j*0.0

            for particle in range(self.particles):
                X_k2=np.array(np.sqrt(s_x[1:self.n])*np.random.normal(size=self.n-1)*np.exp(1j*np.random.rand(self.n-1)*2*np.pi))
                X_k3=np.array(np.random.normal()*np.sqrt(s_x[self.n]))
                X_k1_2=np.append(X_k1,X_k2)
                X_k1_k3=np.append(X_k1_2,X_k3)
                X_k=np.append(X_k1_2,np.conjugate(X_k2[::-1]))
                X_n=np.fft.ifft(X_k)
                y_n=X_n[:self.n]-X_n[0]
                r_n=y_n[:]*np.sqrt(4.*self.D*(self.n)**self.alpha*self.n)
                r_t_allparticles.append(r_n)
                #r_t_allparticles.append(r_x)
            return np.array(r_t_allparticles)


        if self.version=="lowencpp":
            inc1 = ginc1.pyIncrements(self.n,self.particles)
            inc1.generateIncrements1(self.D, self.dt, self.alpha)
            a =inc1.returnIncrements()
            r_t=np.zeros((self.particles,self.n))
            r_t[:,1:]=np.cumsum(a[:,0,1:self.n],axis=1)
            return r_t

        if self.version=="nhoch3":
            r_t_allparticles=[]
            covariance = np.zeros((self.n,self.n))
            A = np.zeros((self.n,self.n))

            for i in range(self.n):
                for j in range(self.n):
                    d = abs(i-j)
                    covariance[i,j] = (abs(d - 1)**self.alpha + (d + 1)**self.alpha - 2*d**self.alpha)*(self.dt**self.alpha*self.D)
            w,v = np.linalg.eig(covariance)
            for i in range(self.n):
                for j in range(self.n):
                    A[i,j] = sum(np.sqrt(w[k])*v[i,k]*v[j,k] for k in range(self.n))
            for particle in range(self.particles):
                x = np.random.randn((self.n))
                eta = np.dot(A,x)
                r_t_allparticles.append(np.cumsum(eta)-eta[0])
            return np.array(r_t_allparticles)

        if self.version=="lowenmodified":
            inc3 = ginc3.pyIncrements(self.n,self.particles)
            inc3.generateIncrements1(self.D, self.dt, self.alpha)
            a =inc3.returnIncrements()
            r_t=np.zeros((self.particles,self.n))
            r_t[:,1:]=np.cumsum(a[:,0,1:self.n],axis=1)
            return r_t

        if self.version=="hartecppverion2":
            inc3 = ginc4.pyIncrements(self.n,self.particles)
            inc3.generateIncrements1(self.D, self.dt, self.alpha)
            a =inc3.returnIncrements()
            r_t=np.zeros((self.particles,self.n))
            r_t[:,1:]=np.cumsum(a[:,0,1:self.n],axis=1)
            return r_t

        if self.version=="lowennew":
            r_t_allparticles=[]
            n=np.array(range(2*self.n))*1.
            r_x1=(self.dt**self.alpha)*(1-(n[:self.n+1]/self.n)**self.alpha)/2
            r_x2=r_x1[::-1]

            r_x=np.append(r_x1,r_x2[1:])
            s_x=np.fft.fft(r_x)
            X_k1=0+0j

            for particle in range(self.particles):
                X_k2=np.array(np.sqrt(s_x[1:self.n]/2)*(np.random.normal(size=self.n-1)+1j*np.random.normal(size=self.n-1)))
                X_k3=np.array(np.random.normal()*np.sqrt(1*s_x[self.n]))
                X_k1_2=np.append(X_k1,X_k2)
                X_k1_k3=np.append(X_k1_2,X_k3)
                X_k=np.append(X_k1_2,np.conjugate(X_k2[::-1]))

                X_n=np.fft.ifft(X_k)
                y_n=X_n[:self.n]-X_n[0]
                r_n=y_n[:]*np.sqrt(4.*self.D*(self.n)**self.alpha*self.n)
                r_t_allparticles.append(r_n)
                #r_t_allparticles.append(r_x)
            return np.array(r_t_allparticles)


