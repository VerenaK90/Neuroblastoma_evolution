This folder contains scripts to analyze and plot the data. 

The file Settings.R should be sourced prior to the analysis. Figures can be reproduced by sourcing the respective script (FigureX.R).

In Settings.R, you need to specify the directories containing the data and the directories, in which the output is to be written. The data should be organized in the same directory with subdirectories named by the tumor IDs. Each tumor subdirectory must contain two folders:
- ACEseq, containing the ACEseq output for this tumor (*comb_pro_extra*)
- SNVs, containing the .vcf file for this tumor. 

In order to fit a model of neutral evolution to a measured VAF, specify the tumor ID in the file Model_neutral_tumor_evolution.R (which contains the functions to simulate a VAF distribution under a neutral model of tumor growth) and execute Neutral_fit.py, which contains the prior probabilities for the parameters and specifies the fiting algorithm.

For a given tumor, you need to specify the purity and the ploidy in the file. Also specify the tumor id and the data directory (for example the Example_data folder). 

The parameters are defined as

n_clonal: # clonal mutations
mu: mutation rate
delta: relative loss rate

The folder contains also scripts for modeling neuroblastoma initiation in two events and estimating the associated parameters (see below for definition) from mutation density data at ECA and MRCA. 

To execute the script, first set your working directory to the folder containing the analysis scripts. Then adjust the file Neuroblastoma_initiation.R. Here, you need to select the growth dynamics of your stem cell population. Set "expansion_decay" for an initially expanding and then contracting population and "expansion_homeostasis" for an initially expanding and subsequently homeostatic population. 

Input data is provided in the file "Input_data.NB.initiation.RData", which is automatically loaded. This file contains the mutation densities at MRCA (NB_origin) and ECA (NB_origin.eca) of tumors with a telomere maintenance mechanism, along with the 95% CI bounds as  estimated in the manuscript. 

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


The script Mutation_density_quantification.R is a pipeline to estimate mutation densities at MRCA and ECA from tumor WG sequencing data. In the first part, mutations are classified as subclonal, non-amplified clonal or amplified clonal. In the second part, these estimates are being used to learn the densities at MRCA and ECA. 

The key output are the vectors ‘mutation.time.mrca’, ‘mutation.time.mrca.lower’,  ‘mutation.time.mrca.upper’,  ‘mutation.time.eca’,  ‘mutation.time.eca.upper’,  ‘mutation.time.eca.lower’, which contain the estimated mutation load per genome copy at MRCA or ECA per tumor (and upper and lower 95% CI according to bootstrapping).
In addition, the following lists are returned 
- ‘gains.at.mrca’: chromosomal gains mapping to MRCA
- ‘gains.uniquely.mapped.to.eca’: chromosomal gains not mapping to MRCA, but to ECA
- ‘gains.at.mrca.conforming.eca’: chromosomal gains mapping to both MRCA and ECA
- ‘gains.not.maping.to.eca.or.mrca’: chromosomal gains neither mapping to MRCA nor to ECA
- ‘gains.at.earliest.time’: chromosomal gains mapping to a time point earlier than ECA


Example data are provided in the folder “Example_data”.
