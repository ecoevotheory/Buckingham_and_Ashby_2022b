This ReadMe file contains information on files associated with the paper: "The evolution of the age of onset of resistance to infectious disease" by Buckingham & Ashby.

All files are provided "as is", without any express or implied warranty. The authors accept no liability for damages of any kind. 

Author: Lydia Buckingham, University of Bath
Contact: ljb74@bath.ac.uk
Date: 13/10/22

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

COMMENTS

All code refers to different combinations of trade-offs as numbered 'versions' of the model, as follows:
version = 1 : Constant fecundity costs paid only when the host is in the resistant stage
version = 2 : Constant mortality costs paid only when the host is in the resistant stage
version = 3 : Trade-off with fecundity throughout the host's lifetime
version = 4 : Trade-off with mortality throughout the host's lifetime
version = 5 : Trade-off with fecundity paid only when the host is susceptible or infected
version = 6 : Trade-off with mortality paid only when the host is susceptible or infected

Figures are generated using the source code below (with parameter values as given in the code):
Figure 1:   "tradeoff_plots_figure.m"
Figure 2:   "constant_costs_singstrat_figure.m"
Figure 3:   "constant_costs_simulation_figure.m"
Figure 4:   "heatmaps_varying_virulence_figure.m"
Figure 5:   "transmissibility_curves_figure.m"
Figure 6:   "branching_simulation_figure.m"
Figure 7:   "tradeoff_shape_effect_figure.m"
Figure S1:  "lifelong_costs_singstrat_figures.m"
Figure S2:  "preresistant_costs_singstrat_figures.m"

MEX files written in C# must be compiled before use, using " mex codename.c ".

Be aware that some of the code relies on numerical approximations and so results may not be exact. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DESCRIPTION OF FILES

Figure 1:
"tradeoff_plots_figure.m"	         - Source code for plotting example trade-off functions									 - written in matlab (R2019b). 


Figure 2:
"constant_costs_singstrat_figure.m"      - Source code for plotting the effect of different parameters on the singular value of zeta				 - written in matlab (R2019b). 


Figure 3:
"constant_costs_simulation_figure.m"     - Source code for plotting evolutionary trajectories for the rate of onset of resistance			         - written in matlab (R2019b)

"simulation_function.m" 		 - Function which runs an evolutionary simulation for a fixed number of timesteps					 - written in matlab (R2019b). 
"Simulation_constantcost_function.c" 	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 1 or 2 costs  - written as a MEX file in C#. 
"Simulation_aTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 3 costs	 - written as a MEX file in C#. 
"Simulation_bTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 4 costs  	 - written as a MEX file in C#. 
"Simulation_earlyaTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 5 costs 	 - written as a MEX file in C#. 
"Simulation_earlybTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 6 costs 	 - written as a MEX file in C#. 


Figure 4:
"heatmaps_varying_virulence_figure.m"    - Source code for drawing heatmaps showing the effect of pathogen virulence on singular strategies			 - written in matlab (R2019b)

"find_singstrats_function_heatmap.m"     - Function which determines singular strategies for given parameter values and determines whether they are CSS's        - written in matlab (R2019b)
"find_singstrats_function2.m"		 - Function which finds singular strategies and determines their evolutionary and convergence stability			 - written in matlab (R2019b)
"classification_simulation_function.m"   - Function which uses simulations to find the stability of singular strategies when it cannot be determined numerically - written in matlab (R2019b)
"simulation_function.m" 		 - Function which runs an evolutionary simulation for a fixed number of timesteps					 - written in matlab (R2019b). 
"Simulation_constantcost_function.c" 	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 1 or 2 costs  - written as a MEX file in C#. 
"Simulation_aTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 3 costs	 - written as a MEX file in C#. 
"Simulation_bTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 4 costs  	 - written as a MEX file in C#. 
"Simulation_earlyaTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 5 costs 	 - written as a MEX file in C#. 
"Simulation_earlybTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 6 costs 	 - written as a MEX file in C#. 


Figure 5:
"transmissibility_curves_figure.m"	 - Source code for plotting the effect of transmissibility on the singular values of rate of onset of resistance	 - written in matlab (R2019b)

"find_singstrats_function.m"		 - Function which determines singular strategies for given parameter values and determines their stability		 - written in matlab (R2019b)
"find_singstrats_function2.m"		 - Function which finds singular strategies and determines their evolutionary and convergence stability			 - written in matlab (R2019b)
"classification_simulation_function.m"   - Function which uses simulations to find the stability of singular strategies when it cannot be determined numerically - written in matlab (R2019b)
"simulation_function.m" 		 - Function which runs an evolutionary simulation for a fixed number of timesteps					 - written in matlab (R2019b). 
"Simulation_constantcost_function.c" 	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 1 or 2 costs  - written as a MEX file in C#. 
"Simulation_aTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 3 costs	 - written as a MEX file in C#. 
"Simulation_bTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 4 costs  	 - written as a MEX file in C#. 
"Simulation_earlyaTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 5 costs 	 - written as a MEX file in C#. 
"Simulation_earlybTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 6 costs 	 - written as a MEX file in C#. 


Figure 6:
"branching_simulation_figure.m"		 - Source code for plotting evolutionary trajectories in cases where branching occurs					 - written in matlab (R2019b)

"simulation_function.m" 		 - Function which runs an evolutionary simulation for a fixed number of timesteps					 - written in matlab (R2019b). 
"Simulation_constantcost_function.c" 	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 1 or 2 costs  - written as a MEX file in C#. 
"Simulation_aTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 3 costs	 - written as a MEX file in C#. 
"Simulation_bTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 4 costs  	 - written as a MEX file in C#. 
"Simulation_earlyaTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 5 costs 	 - written as a MEX file in C#. 
"Simulation_earlybTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 6 costs 	 - written as a MEX file in C#. 


Figure 7:
"tradeoff_shape_effect_figure.m"	 - Source code for plotting the effect of the trade-off shape on singular strategies					 - written in matlab (R2019b)

"find_singstrats_function.m"		 - Function which determines singular strategies for given parameter values and determines their stability		 - written in matlab (R2019b)
"find_singstrats_function2.m"		 - Function which finds singular strategies and determines their evolutionary and convergence stability			 - written in matlab (R2019b)
"classification_simulation_function.m"   - Function which uses simulations to find the stability of singular strategies when it cannot be determined numerically - written in matlab (R2019b)
"simulation_function.m" 		 - Function which runs an evolutionary simulation for a fixed number of timesteps					 - written in matlab (R2019b). 
"Simulation_constantcost_function.c" 	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 1 or 2 costs  - written as a MEX file in C#. 
"Simulation_aTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 3 costs	 - written as a MEX file in C#. 
"Simulation_bTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 4 costs  	 - written as a MEX file in C#. 
"Simulation_earlyaTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 5 costs 	 - written as a MEX file in C#. 
"Simulation_earlybTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 6 costs 	 - written as a MEX file in C#. 


Figure S1:
"lifelong_costs_singstrat_figures.m"     - Source code for plotting the effect of different parameters on singular strategies & their stability (versions 3 & 4) - written in matlab (R2019b)                                                                                    - written in matlab (R2019b). 

"find_singstrats_function.m"		 - Function which determines singular strategies for given parameter values and determines their stability		 - written in matlab (R2019b)
"find_singstrats_function2.m"		 - Function which finds singular strategies and determines their evolutionary and convergence stability			 - written in matlab (R2019b)
"classification_simulation_function.m"   - Function which uses simulations to find the stability of singular strategies when it cannot be determined numerically - written in matlab (R2019b)
"simulation_function.m" 		 - Function which runs an evolutionary simulation for a fixed number of timesteps					 - written in matlab (R2019b). 
"Simulation_constantcost_function.c" 	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 1 or 2 costs  - written as a MEX file in C#. 
"Simulation_aTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 3 costs	 - written as a MEX file in C#. 
"Simulation_bTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 4 costs  	 - written as a MEX file in C#. 
"Simulation_earlyaTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 5 costs 	 - written as a MEX file in C#. 
"Simulation_earlybTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 6 costs 	 - written as a MEX file in C#. 


Figure S2:
"preresistant_costs_singstrat_figures.m" - Source code for plotting the effect of different parameters on singular strategies & their stability (versions 5 & 6) - written in matlab (R2019b)                                                                                    - written in matlab (R2019b). 

"find_singstrats_function.m"		 - Function which determines singular strategies for given parameter values and determines their stability		 - written in matlab (R2019b)
"find_singstrats_function2.m"		 - Function which finds singular strategies and determines their evolutionary and convergence stability			 - written in matlab (R2019b)
"classification_simulation_function.m"   - Function which uses simulations to find the stability of singular strategies when it cannot be determined numerically - written in matlab (R2019b)
"simulation_function.m" 		 - Function which runs an evolutionary simulation for a fixed number of timesteps					 - written in matlab (R2019b). 
"Simulation_constantcost_function.c" 	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 1 or 2 costs  - written as a MEX file in C#. 
"Simulation_aTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 3 costs	 - written as a MEX file in C#. 
"Simulation_bTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 4 costs  	 - written as a MEX file in C#. 
"Simulation_earlyaTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 5 costs 	 - written as a MEX file in C#. 
"Simulation_earlybTradeoff_function.c"	 - Function which runs ecological dynamics for a single evolutionary timestep, for the model using version 6 costs 	 - written as a MEX file in C#. 


Other code included:
"ecological_stability.m"        	 - Code which shows that the endemic equilibrium of the system is linearly stable for a wide range of parameters   	 - written in matlab (R2019b). 
"constant_costs_convergence_stability.m" - Code which shows that the singular strategies in the constant costs case are convergence stable                       - written in matlab (R2019b). 


See code for full description and instructions for use. 