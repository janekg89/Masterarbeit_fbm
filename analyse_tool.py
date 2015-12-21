__author__ = 'janekg89'
import simulation
import numpy as np # module for scientific computing
import matplotlib.pyplot as plt # module for plotting "a la" matlab


class Analyse(simulation.Felix_Method):
    """
        Beinhaltet mehrere Methoden zur Analyse der Trajectorie. Erbt von der Klasse Simulation.Felix_Method .
    """


    def msdanalyt(self):
        return 2*((self.t)**self.alpha)*self.K_alpha


    def msd_ensemble(self):
        r_t_allparticles=simulation.Felix_Method(self.D,self.particles,self.n,self.alpha).compute_trajectory()
        r_t_allparticles=np.array(r_t_allparticles)
        r_t_allparticles_squared=abs(r_t_allparticles)**2
        msd=r_t_allparticles_squared.mean(axis=0)
        std=r_t_allparticles.std(axis=0)
        return  msd, std

    def plotting(self,error=2):
        msd_ensemble,std=Analyse(self.D,self.particles,self.n,self.alpha).msd_ensemble()
        colors=['r','b','g','k','c','w','b','r','g','b','k','c','w','b','r','g','b','k','c','w','bo','ro','go','bo','ko','co','wo','bo']
        #fig=plt.plot(range(msd_ensemble.size), msd_ensemble ,colors[2], label="ensemble msd")
        plt.loglog(self.t,Analyse(self.D,self.particles,self.n,self.alpha).msdanalyt(),":",color=colors[1], label="analytisch D=%f,particles=%d,length=%d,alpha=%f" %(self.D,self.particles,self.n,self.alpha))
        fig=plt.errorbar(range(msd_ensemble.size), msd_ensemble, yerr=error*std,label="Spektrale Methode mit D=%f,particles=%d, length=%d ,alpha=%f, Std=%f" %(self.D,self.particles,self.n,self.alpha,error))
        plt.legend(loc=2)
        plt.xlabel('Steps', fontsize=14)
        plt.ylabel('MSD', fontsize=14)
        return fig



