#import comlex_anal
import complex_sim
#import numpy as np

class Analyse():
    def __init__(self):
        self.particlenumber=20
        self.boxsize=8.0
        self.D=1.0/6.0
        self.R=1.0
        self.alpha=1.0
        self.lambda_plus=1.0
        self.lambda_c=1.0
        self.lambda_minus=1.0
        self.length=2**14


    def analys_protokoll(self):

        """
        self.__init__()

        for self.boxsize in [5,6,7,8.0,9,(std::find(recordingRange.begin(), recordingRange.end(), timeIndex) != recordingRange.end())10,11,12,13,14,15]:
            for index1 in range(300):
                complex1=complex_sim.Sim_Complex(self.D,self.R,self.lambda_plus,self.lambda_minus,self.lambda_c,1.0,self.length,self.boxsize,self.particlenumber,0.05)
                complex1.observables=["all"]
                print index1
                complex1.run()

        self.__init__()

        for self.lambda_plus in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
            for index1 in range(300):
                complex=complex_sim.Sim_Complex(self.D ,self.R,self.lambda_plus,self.lambda_minus,self.lambda_c,1.0,self.length,self.boxsize,self.particlenumber,0.05)
                complex.observables=["all"]
                print index1
                complex.run()

        self.__init__()

        for self.lambda_minus in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
            for index1 in range(300):
                complex=complex_sim.Sim_Complex(self.D ,self.R,self.lambda_plus,self.lambda_minus,self.lambda_c,1.0,self.length,self.boxsize,self.particlenumber,0.05)
                complex.observables=["all"]
                print index1
                complex.run()

        self.__init__()

        for self.lambda_c in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:

            for index1 in range(300):
                complex=complex_sim.Sim_Complex(self.D ,self.R,self.lambda_plus,self.lambda_minus,self.lambda_c,1.0,self.length,self.boxsize,self.particlenumber,0.05)
                complex.observables=["all"]
                print index1
                complex.run()
        self.__init__()

        for self.boxsize in [5,6,7,8,9,10,11,12,13,14,15]:
            for index1 in range(300):
                complex=complex_sim.Sim_Complex(self.D ,self.R,self.lambda_plus,self.lambda_minus,self.lambda_c,0.5,self.length,self.boxsize,self.particlenumber,0.05)
                complex.observables=["all"]

                print index1
                complex.run()
        self.__init__()

        for self.lambda_plus in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:

            for index1 in range(300):
                complex=complex_sim.Sim_Complex(self.D ,self.R,self.lambda_plus,self.lambda_minus,self.lambda_c,0.5,self.length,self.boxsize,self.particlenumber,0.05)
                complex.observables=["all"]

                print index1
                complex.run()
        self.__init__()

        for self.lambda_minus in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:

            for index1 in range(300):
                complex=complex_sim.Sim_Complex(self.D ,self.R,self.lambda_plus,self.lambda_minus,self.lambda_c,0.5,self.length,self.boxsize,self.particlenumber,0.05)
                complex.observables=["all"]

                print index1
                complex.run()
        self.__init__()

        for self.lambda_c in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,]:
            for index1 in range(300):
                complex=complex_sim.Sim_Complex(self.D ,self.R,self.lambda_plus,self.lambda_minus,self.lambda_c,0.5,self.length,self.boxsize,self.particlenumber,0.05)
                complex.observables=["all"]
                print index1
                complex.run()
        self.__init__()
        """

        for self.lambda_minus in [0.0]:
            self.lambda_plus=0.1
            for index1 in range(10):
                complex=complex_sim.Sim_Complex(self.D ,self.R,self.lambda_plus,self.lambda_minus,self.lambda_c,1.0,self.length,8.014,self.particlenumber,0.05)
                complex.observables=["all"]
                print index1
                complex.run()







Analyse().analys_protokoll()
