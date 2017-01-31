## fBm generating algorithms
This directory contains fBm genrating algorithms. simulation.py is containing various fBm generating algorithms. 
The folders davis_harte, lowen, lowen, lowen_modified, naive contain the coresponding C++ algorithms in files "_generatefracincrements.cpp". 
File "analysis_fBm.py" contains an analysis class. 
The ipython notebook "plotting.ipynb" contains performance, accuracy analysis. 
Results from the "plotting.ipynb" are saved in the directory "results".

### Dependencies

Linked libraries:
"fftw3","m","gsl","gslcblas"

### Building
After a modification of C++ code:
The python binding to the C++ code  in "_generatefracincrements.cpp"  in directories "davis_harte, lowen, lowen_modified, naive" can be compilied by:
	
	$ cd folder_of_algorithm/
	$ python setup.py build_ext --inplace 

### Usage

It is recommmended to run the algorithm by file "analysis_fBm.py" as follows:

1. initializing the class Analyse(...) with the desired parameters. The parameter "version" determines the fBm-generating algorithm.
2. The resulting trajecory is safed in self.trajectory.


Example:

import analysis_fBm

lowen=analysis_fBm.Analyse(D=2,particles=5000,length=length+1,alpha=alpha,dt=1.0,version="lowen_modified")

'''
variabe "lowen.trajectory"  contains the trajectory
'''

