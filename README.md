# Neuroblastoma evolution

This repository contains custom code associated with the manuscript "Genetic origin of neuroblastoma in early fetal development" by Verena Körber et al.

To run the analysis, download the R package NBevolution_0.0.0.9000.tar.gz and install it. Thereafter, either refer to the 2 tutorials provided as pdfs in the folder Example to run parts of the analysis on exemplary data or re-run the entire analysis by following the steps explained below. In order to re-run the full analysis you need to download the entire dataset from Mendeley as well as the supplementary tables from the manuscript.

In the folder Analysis_and_plots you can find custom scripts written to analyze the data with the functions provided in the R package "NBevolution" as well as scripts used for plotting. Associated example data (copy number information and SNVs) is provided in the folder Example_data. 

## Contents


- [Software requirements](#software-requirements)
- [Example data](#Example-data)
- [Full dataset](#Full-dataset)
- [Description of folders, data and scripts](#Description-of-folders,-data-and-scripts)
- [Estimate the mutation density at MRCA and ECA for a single example tumor](#Estimate-the-mutation-density-at-MRCA-and-ECA-for-a-single-example-tumor)
- [Analyze the data with population genetics models using uploaded input data](#Analyze-the-data-with-population-genetics-models-using-uploaded-input-data)
- [Reproduce full analysis](#Reproduce-full-analysis)

## Software requirements


R v4.0.0

libraries: 
bedr v1.0.7 
Hmisc v4.4-0
mixtools v1.2.0
scales v1.1.1

python v3.6.1

pyABC v0.9.13

## Example data
The folder Example_data contains SNVs and copy number information for two primary-relapse-sample pairs. Each folder contains two subfolders, ACEseq, storing the copy number information (*comb_pro_extra*) and SNVs, storing the vcf file for somatic SNVs. Both files were generated using DKFZ inhouse-pipelines, which are described under https://github.com/DKFZ-ODCF/SNVCallingWorkflow and under https://github.com/DKFZ-ODCF/ACEseqWorkflow.

## Full dataset
The full dataset (SNVs, Indels, SVs and CNVs for each tumor) can be downloaded from Mendeley (link will be provided upon publication). 

## Description of folders, data and scripts

The package NBevolution contains custom functions used in the analysis. A brief description of each function is provided in the following, a more detailed explanation is given in the help file of each function. The scripts are tailored towards an analysis using SNV/indel output generated by the DKFZ SNV/indel calling pipeline (https://github.com/DKFZ-ODCF/SNVCallingWorkflow) and CNV estimates generated by ACEseq (https://github.com/DKFZ-ODCF/ACEseqWorkflow).

### Extract.info.from.vcf.R 
Extracts information such as the number of reference and variant reads, the variant allele frequency, or the amino acid change from a vcf file as generated by the DKFZ inhouse SNV/indels calling pipeline (https://github.com/DKFZ-ODCF/SNVCallingWorkflow).

### Extract.purity.ploidy.from.ACEseq.R 
Extracts the purity and ploidy estimate suggested by ACEseq (https://github.com/DKFZ-ODCF/ACEseqWorkflow) from the name of the ACEseq *comb_pro_extra* output file. 

### Extract.copy.number.info.per.SSNV.R 
This function assigns to each mutation in a vcf file the copy number state (coverage.ratio (bias corrected coverage ratio), baf (B-allele frequency), genotype (ratio of rounded allele copy numbers), tcn (rounded copy number) and A (rounded minor allele copy number) as suggested by ACEseq. 

### Count.clonal.mutations.R 
This function counts the number of clonal SNVs per copy number state, stratified by mutation multiplicity. The function requires as input the location of the ACEseq CNV estimate as well as of the SNV file. By default, the function analyses clonal SNVs on copy numbers up to four. For each copy number – multiplicity combination, it quantifies the number of SNVs residing at a single or at multiple clonal copies by fitting Binomial distributions to the variant allele frequencies of clonal SNVs and inferring the weights associated with the respective multiplicities. To avoid overestimation of single-copy clonal variants due to putative superposition of the clonal peak with subclonal variants, we estimated the number of single-copy clonal SNVs from the right-hand-side of the peak only. The function returns the ploidy and purity estimate as obtained with ACEseq and two tables giving the number of clonal variants per copy number / multiplicity state per chromosome as well as the corresponding segment lengths.   

### Mutation.time.converter.R
This function converts clonal mutation counts into mutation times. In order to compare mutations on single copies to mutations on multiple copies, the estimates need to be normalized because on gained alleles, single-copy clonal variants only started accumulating after the gain of the chromosomal copies, while on normal alleles, single-copy clonal variants were accumulated all the time. Thus in order to compare the density of single-copy clonal variants between gained and normal segments, the measured number of single-copy clonal variants on gained alleles is extrapolated to “mutation time” by adding the number of multi-copy clonal variants to the number of single-copy clonal variants on each copy. A per-copy mutation density is then obtained by dividing the adjusted number of single-copy clonal variants by the total number of copies.  

### MRCA.ECA.quantification.R
This function estimates the mutation density at the most recent common ancestor (MRCA) of the respective tumor. In addition, it tests whether there is statistical evidence for an early common ancestor (ECA) preceding the MRCA. The function takes two matrices providing the estimated number of clonal mutations per copy number state and the respective segment lengths stratified by chromosome as input and returns the estimated mutation time at ECA and MRCA together with empirical 95% confidence intervals (obtained by bootstrapping the genomic segments). It, moreover, returns the segments mapping to MRCA, ECA along with segments mapping neither to MRCA nor ECA. The mutation times are given as #SNVs per human genome (3.3x109bp).

### Model_NB.initiation.R 
This script contains a set of functions to fit a population genetics model of neuroblatoma initiation to the measured densities at ECA and MRCA. 

## Estimate the mutation density at MRCA and ECA for a single example tumor
Please refer to the tutorial https://github.com/hoefer-lab/Neuroblastoma_evolution/Example/Example_code_mutation_density.pdf and also download the folder https://github.com/hoefer-lab/Neuroblastoma_evolution/RData.

## Analyze the data with population genetics models using uploaded input data
Please refer to the tutorial https://github.com/hoefer-lab/Neuroblastoma_evolution/Example/Dynamics_of_neuroblastoma_initiation.pdf and also download the folder https://github.com/hoefer-lab/Neuroblastoma_evolution/RData.

## Reproduce full analysis
To reproduce the analysis do the following steps:

- Download the example data set from the folder https://github.com/hoefer-lab/Neuroblastoma_evolution/Example_data or the full data set from Mendeley (link will be provided upon publication) and store the files on your system. Each tumor contains four subfolders named “ACEseq” (containing the copy number estimates), “SNVs” (containing the vcf-file for somatic snvs), "Indels" (containing the vcf-file for somatic indels) and "SVs" (containing the vcf files for somatic SVs). On Mendeley, we also stored modeling results and signatures. In addition, you will need the information from Supplementary Table 1 (containing meta data of the samples). 
- Install the package https://github.com/hoefer-lab/Neuroblastoma_evolution/NBevolution_0.0.0.9000.tar.gz and download the folders https://github.com/hoefer-lab/Neuroblastoma_evolution/Analysis_and_plots and github.com/hoefer-lab/Neuroblastoma_evolution/Meta_data.
- Open the file https://github.com/hoefer-lab/Neuroblastoma_evolution/Analysis_and_plots/Settings.R and perform the following steps. 
-	Install the libraries as given in the top section of the file Settings.R. 
-	Moreover, set the directories correctly (second section). 
-	Define the directories under which the data of the discovery set (data.directory.discovery) and the validation set (data.directory.validation) are stored. These directories contain one folder for each tumor of the set (e.g., Example_tumor_1).
-	The subdirectories (e.g., snvs.directory, indel.directory, cnv.directory, etc.) name the subdirectories in which the vcf-files and the output from ACEseq are stored for each tumor. In the example data set, snv.directory=”SNVs/” and cnv.directory=”ACEseq/”. 
-	 The variable "custom.script.directory" stores the path to the folder containing the custom scripts as in the folder “Analysis_and_plots”. 
- The variable "output.directory" stores the name of the directory in which the plots generated by the analysis will be stored (subdirectories will be created within this directory). 
- The variable	"Meta.data" stores the path to the folder containing meta data such as Supplementary Table 1 or gencode_v19_gene_pos.txt, giving the position of the genes in hg19.
-	The variable "rdata.directory" stores the path to the folder where output of the analysis is to be stored as .rdata objects. 
-	Once the previous steps are done, source the Settings.R file.

### Steps to reproduce full analysis:
Source the files in the following order
- Figure_1.R, which sources the following files
-	Gains_and_losses.R which characterizes gains and losses in the individual tumors
-	Driver_mutations.R which summarizes known drivers among the tumors
-	Oncoprint_Fig1.R which generates the oncoprint
-	Mutational_signatures.R which calls and summarizes mutational signatures in the data
- Figure_2.R which sources / uses output generated by the following scripts
- Adjust_purity.R which adjusts the purity in tumors with poor estimates and else takes the output from ACEseq
- Mutation_density_quantification.R which computes the mutation density at each chromosomal segment and at ECA and MRCA of each tumor
- Figure_3.R which sources / uses output generated by the following scripts
- Oncoprint.R which summarizes the mutational landscape inthe discovery cohort
- Driver_mutations.R which summarizes known drivers among the tumors
- Mutation_density_quantification.R which computes the mutation density at each chromosomal segment and at ECA and MRCA of each tumor
- Figure_4.R which sources / uses output generated by the following scripts
- Mutation_density_quantification.R which computes the mutation density at each chromosomal segment and at ECA and MRCA of each tumor
- Survival_analysis.R which analyzes the survival statistics in both cohorts
- Oncoprint.R and Oncoprint_validation.R which summarizes the mutational landscape inthe discovery and validation cohort, respectively
- Figure_5.R wich sources / uses output generated by the following scripts
- Survival_analysis.R which analyzes the survival statistics in both cohorts
- Figure_6.R which sources / uses output generated by the following scripts
- Mobster.R which runs Mobster on the discovery set
- The population genetics model of tumor initiation (run pyABC on Expansion_decay_continuous_evol.R / Expansion_decay_continuous_evol.py / Expansion_homeostasis_continuous_evol.R / Expansion_homeostasis_continuous_evol.py upon generating the input data by sourcing the file Input_data_NB_initiation.R
- The population genetics model of tumor growth (run pyABC on Neutral_fit_pre_clonal_and_clonal.R / Neutral_fit_pre_clonal_and_clonal.py upon generating the input data by sourcing the file Mutation_density_quantification.R)
- Compute_evolutionary_parameters_from_growth_model.R, which summarizes the evolutionary parameters generated by the 2 previous steps.
