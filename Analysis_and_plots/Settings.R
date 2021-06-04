## General settings, such as libraries, samples taken, color codes, etc.
##############################################################################################################################################
## Load libraries that are always needed

library("RColorBrewer")
library(ggplot2); theme_set(theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=8, color="black"),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.75),
                                  axis.ticks = element_line(color = "black", size=0.75),
                                  axis.text = element_text(size=8, color="black")))
library(bedr)
library(openxlsx)
library(pammtools)
library(ComplexHeatmap)
library(ggsignif)
library(dplyr)
library(GenomicRanges)
library(ggbio)
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg19)
library(bedr)
library(ggbeeswarm)
library(ggpubr)
library(survminer)
library(survival)
library(Hmisc)
library(scales)
library(mixtools)
library(reshape2)
library(circlize)
library(gridExtra)
library(MASS)
library(HDInterval)
library(cdata)

##############################################################################################################################################
## set directories
## set to base folder containing the data of the discovery set. 
data.directory.80x <- "./DiscoverySet/"
## set to base folder containing the data of the validation set
data.directory.30x <- "./ValidationSet/"
## subdirectory containing SNV files
snv.directory <- ""
## subdirectory containing indel files
indel.directory <- ""
## subdirectory containing SV files
sv.directory <- ""
## subdirectory containing ACEseq files
cnv.directory <- "ACEseq/"
## directory containing R functions
function.directory <- "./RFunctions/"
## directory containing custom scripts
custom.script.directory <- "./Custom_scripts/"
## output directory to store the plots (base directory, will create subfolders for each figure)
output.directory <- "./Plots/"
## directory containing the meta data (=Supplementary Table 1)
meta.data <- "./MetaData/"
## directory to store .RData files
rdata.directory <- "./RData/"

if(!dir.exists(rdata.directory)){
  dir.create(rdata.directory)
}


##############################################################################################################################################
## source functions

source(paste0(function.directory, "Extract.copy.number.info.per.SSNV.R"))
source(paste0(function.directory, "Extract.purity.ploidy.from.ACEseq.R"))
source(paste0(function.directory, "Extract.info.from.vcf.R"))
source(paste0(function.directory, "Count.clonal.mutations.R"))
source(paste0(function.directory, "MRCA_ECA_quantification.R"))
source(paste0(function.directory, "Mutation.time.converter.R"))
source(paste0(function.directory, "Model_NB_initiation.R"))

##############################################################################################################################################
## Discovery set, meta information and colors

sample.information.80x <- read.xlsx(paste0(meta.data, "Supplementary\ Table1.xlsx"), sheet = 1)
rownames(sample.information.80x) <- sample.information.80x$Tumor_ID

## classify stage 1, 4S, 2 as low-risk, stage 3 as intermediate risk and stage 4 as high-risk
sample.information.80x$ManualScore <- "LR"
sample.information.80x[sample.information.80x$Stage == "3",]$ManualScore <- "IR"
sample.information.80x[sample.information.80x$Stage == "4",]$ManualScore <- "HR"

sample.information.80x$Stage[sample.information.80x$Stage %in% c("2.1", "2.2")] <- "2"

tumors.80x <- sample.information.80x$Tumor_ID

## Stratify by Lesion type (merge primary and metastasis)
primary.tumors.80x <- sample.information.80x$Tumor_ID[sample.information.80x$Location %in% c("Primary", "Metastasis")]
relapse.tumors.80x <- sample.information.80x$Tumor_ID[!sample.information.80x$Location %in% c("Primary", "Metastasis")]

## Stratify by molecular subgroup
telomere.classification.80x <- sample.information.80x$Telomere.maintenance.mechanism
names(telomere.classification.80x) <- sample.information.80x$Tumor_ID

## Stratify by ploidy
diploid.tumors.80x <- rownames(sample.information.80x[sample.information.80x$Ploidy==2,])
triploid.tumors.80x <- rownames(sample.information.80x[sample.information.80x$Ploidy==3,])
tetraploid.tumors.80x <- rownames(sample.information.80x[sample.information.80x$Ploidy==4,])
diploid.tumors.80x <- rownames(sample.information.80x[sample.information.80x$Ploidy==2,])


##############################################################################################################################################
## Validation set, sequenced at 30x

sample.information.30x <-  read.xlsx(paste0(meta.data, "Supplementary\ Table1.xlsx"), sheet = 2)
rownames(sample.information.30x) <- sample.information.30x$Tumor_ID

tumors.30x <- rownames(sample.information.30x)

sample.information.30x$ManualScore <- "LR"
sample.information.30x[sample.information.30x$Stage == "3",]$ManualScore <- "IR"
sample.information.30x[sample.information.30x$Stage == "4",]$ManualScore <- "HR"

## Stratify by ploidy
diploid.tumors.30x <- rownames(sample.information.30x[sample.information.30x$Rounded.ploidy==2,])
triploid.tumors.30x <- rownames(sample.information.30x[sample.information.30x$Rounded.ploidy==3,])
tetraploid.tumors.30x <- rownames(sample.information.30x[sample.information.30x$Rounded.ploidy==4,])
diploid.tumors.30x <- rownames(sample.information.30x[sample.information.30x$Rounded.ploidy==2,])

##############################################################################################################################################
## Color palette: 
## early clonal, late clonal, subclonal

manual.colors <- c("#4FB12B", "#176A02") #brewer.pal(9, "PRGn")[c(7,9)]
names(manual.colors) <- c("Late", "Early")

## Sample type
time.colors <- c(Primary="lightgrey", "Relapse" = "darkblue", "Relapse tumor"="darkblue", "Relapse metastasis" = "orange", Metastasis="firebrick", `Relapse 3`="darkblue")

## Ploidy
ploidy.cols <-brewer.pal(4, "Greys")
names(ploidy.cols) <- c("2", "3", "4", "Other")

## Molecular subgroup
telomere.colors <- c(ALT = "#d7191c",
                     TERT = "#ffffbf",
                     MNA = "#abdda4", 
                     Multiple = "#2b83ba",
                     None = "lightgrey")

clinical.risk.colors <- c(HR = "firebrick", IR = "darkgrey", observation = "lightgrey", LR="lightgrey")

stage.colors <-  brewer.pal(5, "Dark2")
names(stage.colors) <- c("1", "2", "3", "4", "4S")
stage.colors["1"] <- "#1A9F77"
stage.colors["4S"] <- "#A5D062"

##############################################################################################################################################
## Driver genes -- according to Ackermann et al., Science 2018 + MYCN, TERT, ATRX

driver.genes <- c("ALK", "ATM", "BRAF", "CCND1", "CDK4", "CDKN2A", "CREBBP", "FGFR1", "LIN28B", "MDM2", "MDM4", "NF1", "PTPN11", 
                  "HRAS", "KRAS", "NRAS", "RRAS", "TP53", "MYCN", "TERT", "ATRX")


## genetic positions 
gene.positions <- read.delim(paste0(meta.data, "gencode_v19_gene_pos.txt"), header=F, row.names = 1)
driver.genes.with.genetic.pos <- data.frame(GENE=unique(driver.genes), gene.positions[unique(driver.genes),])
colnames(driver.genes.with.genetic.pos) <- c("GENE", "CHROM", "START", "END")

### known non-neutral mutations
known.driver.positions <- read.delim(paste0(meta.data, "Driver_mutations_w_position_DoCM_29_Aug_2020.tsv"),
                                     stringsAsFactors = F, sep="\t")



