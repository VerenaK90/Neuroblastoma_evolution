This folder contains scripts to simulate an expected VAF distribution under neutral evolution. This file Model_neutral_tumor_evolution.R contains functions to simulate the distribution at a give parameter set. 

The file Neutral_fit_pre_clonal_and_clonal.R contains data preprocessing and the input functions for pyABC. For a given tumor, you need to specify the purity and the ploidy in the file. Also specify the tumor id and the data directory (for example the Example_data folder). 

Parameter estimation is done by executing Neutral_fit.py, where also the prior probabilities for the parameters can be adjusted. The parameters are defined as

n_clonal: # clonal mutations
mu: mutation rate
delta: relative loss rate