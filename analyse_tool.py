__author__ = 'janekg89'
import simulation
import numpy as np # module for scientific computing
import matplotlib.pyplot as plt # module for plotting "a la" matlab
from scipy.stats import norm
from scipy.stats import moment


class Analyse(simulation.Felix_Method):
    """
        Beinhaltet mehrere Methoden zur Analyse der Trajectorie. Erbt von der Klasse Simulation.Felix_Method .
    """
    def __init__(self, D , particles,length, alpha,dt):
        simulation.Felix_Method.__init__(self,D=D, particles=particles,length=length,alpha=alpha,dt=dt)
        self.trajectory=simulation.Felix_Method(D,particles,length,alpha,dt=dt).compute_trajectory()
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
        :param histpoints: amount of timepoints
        :return: Distribution of particles at time t
        """
        r_t_allparticles=self.trajectory
        r_tt=r_t_allparticles.T
        hist, hist_p = np.histogram(r_tt[t], histpoints, normed=True)
        hist_p=hist_p[1:]
        return hist_p,hist

    def distribution_abs(self,t=0,histpoints=70):
        """
        :param t: point in time
        :param histpoints: amount of timepoints
        :return: Distribution of particles at time t
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
        #r,dist=Analyse(self.D,self.particles,self.n,self.alpha).distribution_abs(t,histpoints)
        res_r=hist_p*(2*self.D*t)**(-self.alpha/2.)
        res_dis=hist*res_r*np.sqrt(2*self.D*t**self.alpha)

        return res_r, res_dis



    def analytical_distribution_of_particles(self,t, r_dis=50):

        r=np.arange(-r_dis, r_dis, 0.01)
        #r=abs(r)
        distrib=np.exp(-(r**2)/(2*2*self.D*t**self.alpha))/np.sqrt(np.pi*2*2*self.D*t**self.alpha)
        distrib1=norm.pdf(r,0,np.sqrt(2*self.D*t**self.alpha))
        return r, distrib, distrib1

    def rescaled_analytical_distribution(self, t, r_dis= 30):
        """

        :param t: point in time (return should be independed of t, but i left possiblitiy to also demonstrate that this is true)
        :param r_dis: the distance which want to be sampled. (not recommended to change)
        :return: analytical rescaled function ( fill in math)
        """
        r=np.arange(0, r_dis, 0.01)
        res_r=r*(2*self.D*t)**(-self.alpha/2.)
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



    def nongaussian_parameter(self):
        moment2poten2=(moment(self.trajectory,moment=2,axis=0))**2
        moment4=moment(self.trajectory,moment=4,axis=0)
        nongaussianparamter=(1/3.)*moment4/moment2poten2-1
        return  nongaussianparamter

    def invert_time(self):
        """
        Dreht die Richtung der Zeit

        :return:
        """
        traject=np.fliplr(self.trajectory)
        traject=np.subtract(traject.T,traject[:,0].T)
        self.trajectory=traject.T
    def plot_freq(self):
            plt.plot(self.frq)
            plt.show()

    def plotting(self, msdtype="ensemble", particlemsdtime=0,error=0, showlegend=None,scale="loglog"):
        """
        :param error: Number of standard deviations from mean, which is shown in the figure


        :return: A figure with plotting of the Ensemble MSD
        """
        if msdtype=="ensemble":
            msd,std=self.msd_ensemble()
        if msdtype=="time":
            msd,std=self.msd_time(particlemsdtime)


        colors=['r','b','g','k','c','w','b','r','g','b','k','c','w','b','r','g','b','k','c','w','bo','ro','go','bo','ko','co','wo','bo']
        #fig=plt.plot(range(msd_ensemble.size), msd_ensemble ,colors[2], label="ensemble msd")
        if scale == "lin":
            plt.plot(self.t*self.dt,self.msdanalyt(),":",color=colors[1], label="analytisch D=%f,particles=%d,length=%d,alpha=%f" %(self.D,self.particles,self.n,self.alpha))
        if scale == "loglog":
            plt.loglog(self.t*self.dt,self.msdanalyt(),":",color=colors[1], label="analytisch D=%f,particles=%d,length=%d,alpha=%f" %(self.D,self.particles,self.n,self.alpha))
        fig=plt.errorbar(self.t*self.dt, msd, yerr=error*std,label="Spektrale Methode mit D=%f,particles=%d, length=%d ,alpha=%f, Std=%f" %(self.D,self.particles,self.n,self.alpha,error))
        if showlegend is not None:
            plt.legend(loc=2)
        plt.xlabel('Steps', fontsize=14)
        plt.ylabel('MSD', fontsize=14)
        return fig
