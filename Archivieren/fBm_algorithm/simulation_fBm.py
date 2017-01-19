# -*- coding: utf-8 -*-
__author__ = 'janek'
import numpy as np # module for scientific computing
import naive.genereatefracincrements as ginc
import lowen_modified.genereatefracincrements as ginc1
import lowen.genereatefracincrements as ginc3
import davis_harte.genereatefracincrements as ginc4

class Simulation_fbm():
    def __init__(self, D, particles,length,alpha,dt,x=2, version="hartecppverion2"):
        """

        :param D: diffusion coeffcient
        :param particles:
        :param length:
        :param alpha:
        :param dt:
        :param x: (only important
        :param version:
        :return:
        """
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

    def z_velocity(self):
        """
        :param D: Diffusionscoefficient
        :param alpha: anomalous diffusion exponent
        :return z_velocity: Velocity-Autocorrelation function (F. Hoefling and T. Franosch: Anomalous Transport in the Crowded World of Biological Cells, Reports on Progress in Physics 76, 046602 (2013))
        """
        z=((((+1j*self.frq))**(1-self.alpha))*self.K_alpha*np.math.gamma(1+self.alpha))
        return z
    def z_fgn(self):
        """
        :param D: Diffusionscoefficient
        :param alpha: anomalous diffusion exponent
        :return z_velocity: fractional Gaussian noise-Autocorrelation function
        """
        z=((((+1j*self.frq1))**(1-self.alpha))*self.K_alpha*np.math.gamma(1+self.alpha))
        return z

    def z_msd(self):
        """
        :return: velocity autocorrelation function calcualted from meansquare displacement
        """
        MSD=(self.t*self.dt)**(self.alpha)*self.K_alpha*2
        # second derivative of fBm MSD
        d2MSD=self.alpha*(self.alpha-1)*(np.array(range(self.ki))*self.dt)**(self.alpha-2.)*self.K_alpha
        z_omega=np.fft.fft(d2MSD[:])
        return z_omega
    def z_alterantive_2(self):
        cov=[]
        for i in self.t[1:-2]:
            cov.append(((i-1)**self.alpha-2*i**self.alpha+(i+1)**self.alpha)/(2))
        cov=np.array(cov)
        return cov



    def compute_trajectory(self):
        """
        :param D: Diffusionskoeffizient
        :param particles: number of particle
        :param length: length of trajectory
        :param alpha: anomoalous diffusion exponent
        :return: ( particles x length ) Array of trajectories of all particles.
        """
        r_t_allparticles=[]
        if self.version == "naive_cpp":
            """
            A naive alogirthm in C++
            """
            inc = ginc.pyIncrements(self.n,self.particles)
            inc.generateIncrements(self.D, self.dt, self.alpha)
            a =inc.returnIncrements()
            r_t=np.zeros((self.particles,self.n))
            r_t[:,1:]=np.cumsum(a[:,0,0:self.n-1],axis=1)
            return r_t
        elif self.version == "naive_python":
            """
            A naive alogirthm in Python
            """
            sqrt2zreal=np.sqrt(self.z_velocity().real*2.)

            for particle in range(self.particles):
                v_t=(np.random.normal(0,np.sqrt(self.dt),size=self.ki))
                v_t=np.array(v_t)
                v_frq=np.fft.fft(v_t)
                v_ano_frq= sqrt2zreal*v_frq
                v_ano_frq[0]=np.random.normal(0,np.sqrt(2.*self.K_alpha*(self.ki*self.dt)**self.alpha))
                v_ano_t=np.fft.ifft(v_ano_frq)
                r_t=np.zeros(self.n)
                r_t[1:]=np.cumsum(v_ano_t[:self.n].real)[:self.n-1]
                r_t_allparticles.append(r_t)
            return np.array(r_t_allparticles)


        elif self.version == "python_notrick":
            """
            For a demonstration in the masterthesis.
            """
            sqrt2zreal=np.sqrt(self.z_fgn().real*2.)
            for particle in range(self.particles):
                v_t=(np.random.normal(0,np.sqrt(self.dt),size=self.n))
                v_t=np.array(v_t)
                v_frq=np.fft.fft(v_t)
                v_ano_frq= sqrt2zreal*v_frq
                v_ano_t=np.fft.ifft(v_ano_frq)
                r_t=np.zeros(self.n)
                r_t[1:]=np.cumsum(v_ano_t[:self.n].real)[:self.n-1]
                r_t_allparticles.append(r_t)
            return np.array(r_t_allparticles)

        elif self.version == "python_trick":
            """
            For a demonstration in the masterthesis.
            """
            sqrt2zreal=np.sqrt(self.z_fgn().real*2.)
            for particle in range(self.particles):
                #r = np.random.RandomState(1234)
                v_t=(np.random.normal(0,np.sqrt(self.dt),size=self.n))
                v_t=np.array(v_t)
                v_frq=np.fft.fft(v_t)
                v_ano_frq= sqrt2zreal*v_frq
                v_ano_frq[0]=np.random.normal(0,np.sqrt(2.*self.K_alpha*(self.ki*self.dt)**self.alpha))
                v_ano_t=np.fft.ifft(v_ano_frq)
                r_t=np.zeros(self.n)
                r_t[1:]=np.cumsum(v_ano_t[:self.n].real)[:self.n-1]
                r_t_allparticles.append(r_t)
            return np.array(r_t_allparticles)

        elif self.version=="lowen":
            """
            Lowen algorithm in python


            S. B. Lowen: Efficent generation of fractional brownian motion for simulation
            of infrared focal-plane array calibration drift, Methodology And Computing
            In Applied Probability 1, 445 (1999)
            """
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
            return np.array(r_t_allparticles)


        if self.version=="lowen_modified":
            """
            Modified lowen algorithm in C++


            S. B. Lowen: Efficent generation of fractional brownian motion for simulation
            of infrared focal-plane array calibration drift, Methodology And Computing
            In Applied Probability 1, 445 (1999)
            """
            inc1 = ginc1.pyIncrements(self.n,self.particles)
            inc1.generateIncrements1(self.D, self.dt, self.alpha)
            a =inc1.returnIncrements()
            r_t=np.zeros((self.particles,self.n))
            r_t[:,1:]=np.cumsum(a[:,0,1:self.n],axis=1)
            return r_t

        if self.version=="pseudo-cholesky":

            """
            Cholesky algorithm
            """

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

        if self.version=="lowencpp":
            """
            Lowen algorithm in C++


            S. B. Lowen: Efficent generation of fractional brownian motion for simulation
            of infrared focal-plane array calibration drift, Methodology And Computing
            In Applied Probability 1, 445 (1999)
            """
            inc3 = ginc3.pyIncrements(self.n,self.particles)
            inc3.generateIncrements1(self.D, self.dt, self.alpha)
            a =inc3.returnIncrements()
            r_t=np.zeros((self.particles,self.n))
            r_t[:,1:]=np.cumsum(a[:,0,1:self.n],axis=1)
            return r_t

        if self.version=="davis-harte":

            """
            Davis-Harte  algorithm in C++


            R. B. DAVIES and D. S. HARTE: Tests for hurst effect, Biometrika 74, 95 (1987)
            """

            inc3 = ginc4.pyIncrements(self.n,self.particles)
            inc3.generateIncrements1(self.D, self.dt, self.alpha)
            a =inc3.returnIncrements()
            r_t=np.zeros((self.particles,self.n))
            r_t[:,1:]=np.cumsum(a[:,0,1:self.n],axis=1)
            return r_t


