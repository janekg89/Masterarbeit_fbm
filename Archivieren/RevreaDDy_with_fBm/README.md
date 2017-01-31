## fBm generating algorithms
This directory contains my implementation of an enzymatic reaction with fractional Brownian motion. The directory "RevreaDDy" contains the modified RevReaDDy simulation tool. 
The README.md in directory "RevreaDDy" explaines the simualtion tool. Fractional Brownian motion can be activated by: s = revreaddy.Sim("FractionalDiffusion") instead of s=rdy.Sim().
The python script "python_model" sets up the enzymatic simulation model. There, the path to the simulated data is defiened (default:"./revreaddy_data"). 
The "simulation_protokoll" simulates 500 trajectories of interessing results for the analysis part. Simulations have to be performed before running "analyse_protokoll". 
It creates plots and saves them into the directory "results". Therefore python script "analysis.py", which  has all analysis functions defiened. 

