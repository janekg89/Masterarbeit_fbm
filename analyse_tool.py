__author__ = 'janekg89'
import simulation
import numpy as np # module for scientific computing
import matplotlib.pyplot as plt # module for plotting "a la" matlab
from scipy.stats import norm


class Analyse(simulation.Felix_Method):
    """
        Beinhaltet mehrere Methoden zur Analyse der Trajectorie. Erbt von der Klasse Simulation.Felix_Method .
    """
    def __init__(self, D , particles,length, alpha):
        simulation.Felix_Method.__init__(self,D=D, particles=particles,length=length,alpha=alpha)
        self.trajectory=simulation.Felix_Method(D,particles,length,alpha).compute_trajectory()
    def msdanalyt(self):
        """
        :return: returns the analytical result for MSD
        .. math::

            \\delta r^{2}(t)=2 D t^{ \\alpha}
        """
        return 2*((self.t)**self.alpha)*self.K_alpha

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
        :param histpoints: amount of timepoints
        :return: Distribution of particles at time t
        """
        r_t_allparticles=self.trajectory
        r_tt=r_t_allparticles.T
        hist, hist_p = np.histogram(r_tt[t], histpoints, normed=True)
        hist_p=hist_p[1:]
        return hist_p,hist


    def rescaled_function(self,t=0,histpoints=70):

        r,dist=Analyse(self.D,self.particles,self.n,self.alpha).distribution(t,histpoints)
        res_r=abs(r)*(2*self.D*t)**(-self.alpha/2.)
        res_dis=dist*res_r*np.sqrt(2*self.D*t**self.alpha)

        return res_r, res_dis



    def analytical_distribution_of_particles(self,t, r_dis=50):
        r=np.arange(-r_dis, r_dis, 0.01)
        #r=abs(r)
        distrib=np.exp(-(r**2)/(2*2*self.D*t**self.alpha))/np.sqrt(np.pi*2*2*self.D*t**self.alpha)
        distrib1=norm.pdf(r,0,np.sqrt(2*self.D*t**self.alpha))
        return r, distrib, distrib1

    def rescaled_analytical_distribution(self, t, r_dis= 50):
        r=np.arange(-r_dis, r_dis, 0.01)
        res_r=abs(r)*(2*self.D*t)**(-self.alpha/2.)
        res_dis=(res_r*np.exp(-res_r**2./2.))/np.sqrt(2.*np.pi)

        #res_r, res_dis = Analyse(self.D,self.particles,self.n,self.alpha).analytical_distribution_of_particles(t,r_dis=r_res)
        #res_r = np.arange(-50, 50, 0.01)
        #res_dis= norm.pdf(res_r,0,np.sqrt(2*self.D))
        return res_r, res_dis

    def msd_time(self,i=0):
        """

        :param i: selected Particle trajectory of the which the time averaged MSD  should be calculated. Cannot be bigger than total partical number

        :return: The time averaged mean square displacement ans the sample standard deviation
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
            std[j]=np.array(msdink).std(axis=0)
        return np.array(msd),np.array(std)

    def plotting(self, msdtype="ensemble", particlemsdtime=0,error=0, showlegend=None):
        """
        :param error: Number of standard deviations from mean, which is shown in the figure


        :return: A figure with plotting of the Ensemble MSD
        """
        if msdtype=="ensemble":
            msd,std=Analyse(self.D,self.particles,self.n,self.alpha).msd_ensemble()
        if msdtype=="time":
            msd,std=Analyse(self.D,self.particles,self.n,self.alpha).msd_time(particlemsdtime)


        colors=['r','b','g','k','c','w','b','r','g','b','k','c','w','b','r','g','b','k','c','w','bo','ro','go','bo','ko','co','wo','bo']
        #fig=plt.plot(range(msd_ensemble.size), msd_ensemble ,colors[2], label="ensemble msd")
        plt.loglog(self.t,Analyse(self.D,self.particles,self.n,self.alpha).msdanalyt(),":",color=colors[1], label="analytisch D=%f,particles=%d,length=%d,alpha=%f" %(self.D,self.particles,self.n,self.alpha))
        fig=plt.errorbar(range(msd.size), msd, yerr=error*std,label="Spektrale Methode mit D=%f,particles=%d, length=%d ,alpha=%f, Std=%f" %(self.D,self.particles,self.n,self.alpha,error))
        if showlegend is not None:
            plt.legend(loc=2)
        plt.xlabel('Steps', fontsize=14)
        plt.ylabel('MSD', fontsize=14)
        return fig