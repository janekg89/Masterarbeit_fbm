__author__ = 'janekg89'
import simulation
import numpy as np # module for scientific computing
import matplotlib.pyplot as plt # module for plotting "a la" matlab
import scipy.stats.mstats

class Analyse(simulation.Felix_Method):
    """
        Beinhaltet mehrere Methoden zur Analyse der Trajectorie. Erbt von der Klasse Simulation.Felix_Method .
    """


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

            \\delta r^{2}(t=\\frac{\\sum_{j=0}^N {\\eta_i}}{N}

            N= number of particles

            i= cumulant

        """
        r_t_allparticles=simulation.Felix_Method(self.D,self.particles,self.n,self.alpha).compute_trajectory()
        r_t_allparticles=np.array(r_t_allparticles)
        r_t_allparticles_squared=abs(r_t_allparticles)**2
        msd=r_t_allparticles_squared.mean(axis=0)
        std=r_t_allparticles.std(axis=0)
        return  msd, std

    def msd_time(self,i=0):
        sinfgletrajectory=simulation.Felix_Method(self.D,self.particles,self.n,self.alpha).compute_trajectory()[i]
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

    def plotting(self,error=2, msdtype="ensemble"):
        """
        :param error: Number of standard deviations from mean, which is shown in the figure


        :return: A figure with plotting of the Ensemble MSD
        """
        if msdtype=="ensemble":
            msd,std=Analyse(self.D,self.particles,self.n,self.alpha).msd_ensemble()
        if msdtype=="time":
            msd,std=Analyse(self.D,self.particles,self.n,self.alpha).msd_time()


        colors=['r','b','g','k','c','w','b','r','g','b','k','c','w','b','r','g','b','k','c','w','bo','ro','go','bo','ko','co','wo','bo']
        #fig=plt.plot(range(msd_ensemble.size), msd_ensemble ,colors[2], label="ensemble msd")
        plt.plot(self.t,Analyse(self.D,self.particles,self.n,self.alpha).msdanalyt(),":",color=colors[1], label="analytisch D=%f,particles=%d,length=%d,alpha=%f" %(self.D,self.particles,self.n,self.alpha))
        fig=plt.errorbar(range(msd.size), msd, yerr=error*std,label="Spektrale Methode mit D=%f,particles=%d, length=%d ,alpha=%f, Std=%f" %(self.D,self.particles,self.n,self.alpha,error))
        plt.legend(loc=2)
        plt.xlabel('Steps', fontsize=14)
        plt.ylabel('MSD', fontsize=14)
        return fig