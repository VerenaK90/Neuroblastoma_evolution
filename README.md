# Neuroblastoma_evolution

This repository contains custom code associated with the manuscript "Evolutionary dynamics of neuroblastoma" by Verena Körber et al.

Folder Mutation_density contains scripts for estimating the mutation densities at MRCA and ECA. Associated example data (copy number information and SNVs) is provided in the folder Example_data.

Folder Modeling_tumor_initiation contains scripts to fit a model of neuroblastoma initiation to estimated mutation densities at MRCA and ECA. 

Folder Modeling_tumor_progression contains scripts to fit a model of neutral mutation accumulation to a measured VAF-distribution. 

Software requirements:

R v4.0.0

libraries: 
bedr v1.0.7 
Hmisc v4.4-0
mixtools v1.2.0
scales v1.1.1

python v3.6.1

pyABC v0.9.13
