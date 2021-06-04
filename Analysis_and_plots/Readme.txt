This folder contains scripts to analyze and plot the data. 

The file Settings.R should be sourced prior to the analysis. Figures can be reproduced by sourcing the respective script (FigureX.R).

In order to fit a model of neutral evolution to a measured VAF, specify the tumor ID in the file Model_neutral_tumor_evolution.R (which contains the functions to simulate a VAF distribution under a neutral model of tumor growth) and execute Neutral_fit.py, which contains the prior probabilities for the parameters and specifies the fiting algorithm.

For a given tumor, you need to specify the purity and the ploidy in the file. Also specify the tumor id and the data directory (for example the Example_data folder). 

The parameters are defined as

n_clonal: # clonal mutations
mu: mutation rate
delta: relative loss rate
