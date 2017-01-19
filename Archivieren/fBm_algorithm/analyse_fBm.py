__author__ = 'janekg89'

import simulation_fBm
import numpy as np # module for scientific computing
from scipy.stats import norm
from scipy.stats import moment


class Analyse(simulation_fBm.Simulation_fbm):
    """
    contains several functions to analyse fbm generating trajectories
    """
    def __init__(self, D , particles,length, alpha,dt,x=2, version="python"):
        simulation_fBm.Simulation_fbm.__init__(self,D=D, particles=particles,length=length,alpha=alpha,dt=dt,x=x,version=version)
        self.trajectory=simulation_fBm.Simulation_fbm(D,particles,length,alpha,dt=dt,x=x,version=version).compute_trajectory()

    def msdanalyt(self):
        """
        :return: returns the analytical result for MSD
        .. math::

            \\delta r^{2}(t)=2 D t^{ \\alpha}
        """
        return 2*((self.t*self.dt)**self.alpha)*self.K_alpha

    def msd_ensemble(self):
        """

        :return: returns the ensemble average of the simulated trajectories of all particles.
        .. math::

            \\delta r^{2}(t)=\\frac{\\sum_{j=0}^N {\\eta_i}}{N}

            N= number of particles

            i= cumulant

        """
        r_t_allparticles=self.trajectory
        r_t_allparticles_squared=abs(r_t_allparticles)**2
        msd=r_t_allparticles_squared.mean(axis=0)
        std=r_t_allparticles.std(axis=0)
        return  msd, std

    def distribution(self,t=0,histpoints=70):
        """
        :param t: point in time
        :param histpoints: number points in histogram
        :return: distribution of particles at time t
        """
        r_t_allparticles=self.trajectory
        r_tt=r_t_allparticles.T
        hist, hist_p = np.histogram(r_tt[t], histpoints, normed=True)
        hist_p=hist_p[1:]
        return hist_p*self.dt,hist

    def distribution_abs(self,t=0,histpoints=70):
        """
        :param t: point in time
        :param histpoints: number points in histogram
        :return:  distribution of particles at time t
        """
        r_t_allparticles=abs(self.trajectory)
        r_tt=r_t_allparticles.T
        hist, hist_p = np.histogram(r_tt[t], histpoints, normed=True)
        hist_p=hist_p[1:]
        return hist_p,hist

    def rescaled_function(self,t=0,histpoints=70):
        r_t_allparticles=abs(self.trajectory)
        r_tt=r_t_allparticles.T
        hist, hist_p = np.histogram(r_tt[t], histpoints, normed=True)
        hist_p=hist_p[1:]
        res_r=hist_p*(2*self.D*(t*self.dt)**(self.alpha))**(-1./2.)
        res_dis=hist*hist_p/2

        return res_r, res_dis



    def analytical_distribution_of_particles(self,t, r_dis=50):

        r=np.arange(-r_dis, r_dis, 0.01)
        distrib=np.exp(-(r**2)/(2*2*self.D*(t*self.dt)**self.alpha))/np.sqrt(np.pi*2*2*self.D*(t*self.dt)**self.alpha)
        distrib1=norm.pdf(r,0,np.sqrt(2*self.D*(t*self.dt)**self.alpha))
        return r, distrib, distrib1

    def rescaled_analytical_distribution(self, t, r_dis= 4):
        """
        :param t: point in time (return should be independed of t, but i left possiblitiy to also demonstrate that this is true)
        :param r_dis: the distance which want to be sampled. (not recommended to change)
        :return: analytical rescaled function
        """
        r=np.arange(0.07, r_dis, 0.01)
        res_dis=r*(np.exp(-r**2./2.))/(np.sqrt(2.*np.pi))

        return r, res_dis

    def msd_time(self,i=0):
        """
        :param i: selected particle trajectory of the which the time averaged MSD should be calculated. Cannot be bigger than total partical number
        :return: The time averaged mean square displacement and the sample standard deviation
        """
        sinfgletrajectory=self.trajectory[i]
        totalsize=len(sinfgletrajectory)
        msd=[None]*len(sinfgletrajectory)
        std=[None]*len(sinfgletrajectory)
        for j in range(0,totalsize):
            msdink=[]
            for i in range(0,totalsize-j):
                msdink.append(((sinfgletrajectory[i]-sinfgletrajectory[j+i])**2))
            msd[j]=np.array(msdink).mean(axis=0)
            std[j]=np.array(msdink).std(axis=0)/np.sqrt(totalsize-j)
        return np.array(msd),np.array(std)



    def nongaussian_parameter(self):
        moment2poten2=(moment(self.trajectory,moment=2,axis=0))**2
        moment4=moment(self.trajectory,moment=4,axis=0)
        nongaussianparamter=(1/3.)*moment4/moment2poten2-1
        return  nongaussianparamter

    def invert_time(self):
        """
        reverse time

        :return:
        """
        traject=np.fliplr(self.trajectory)
        traject=np.subtract(traject.T,traject[:,0].T)
        self.trajectory=traject.T


