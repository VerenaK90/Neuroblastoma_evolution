This folder contains scripts for modeling neuroblastoma initiation in two events and estimating the associated parameters (see below for definition) from mutation density data at ECA and MRCA. 

To execute the script, first set your working directory to the folder containing the analysis scripts. Then adjust the file Neuroblastoma_initiation.R. Here, you need to select the growth dynamics of your stem cell population. Set "expansion_decay" for an initially expanding and then contracting population and "expansion_homeostasis" for an initially expanding and subsequently homeostatic population. 

Input data is provided in the file "Input_data.RData", which is automatically loaded in. This file contains the mutation densities at MRCA (NB_origin) and ECA (NB_origin.eca) of high-risk tumors and at the MRCA of low-risk tumors (NB_origin.mrca.lr), along with the 95% CI bounds as  estimated in the manuscript. 

For parameter estimation execute the script Neuroblastoma_initiation.py (You need to install pyABC, https://pyabc.readthedocs.io/en/latest/, to do this). Prior probabilities can also be adjusted in that script. The parameters encode the following:

delta1: relative loss rate during expansion
N: maximum stem cell count (log10)
delta2: relative loss rate during decay
psurv: survival probability associated with the second oncogenic event
muD1: mutation rate per cell division associated with the first oncogenic event (log10)
muD2: mutation rate per division associated with the second oncogenic event (log10)
mu: mutation rate per cell division (neutral + driver mutations)
r: selective advantage associated with the fist oncogenic event.

If you run "expansion_homeostasis", remove the prior for delta2 from Neuroblastoma_initiation.py.
