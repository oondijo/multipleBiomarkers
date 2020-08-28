# multipleBiomarkers
Simulation code (R scripts) for simulation study comparing treatment allocation strategies for umbrella trials in the presence of multiple biomarkers

Files contained in this repository can be used to reproduce the numerical results reported in the paper entitled
# LO Ondijo, MJ Grayling, H Zheng, JMS Wason. (2020) Treatment allocation strategies for umbrella trials in the presence of multiple biomarkers: A comparison of methods.
*A submitted paper*.

The files "BAR_code_multipleBiomarkers.R", "inits_K=4.txt" and "model_K=4_betterthancontrol" are used to implement the BAR treatment allocation procedure.
The files "ER_RFAC_CR_Hier_functions.R", "ER_RFAC_CR_Hier_Operational_xtics.R", and "ER_RFAC_CR_Hier_Operational_xtics_SuperFn.R" are used to evaluate the 
operational characteristics of 4 approaches namely: Hierarchy, Constrained randomization, Equal randomization and Randomization with fixed allocation probability to control
The file "simulation_scenarios.R" contains the 12 simulation scenarios evaluated in the paper.
