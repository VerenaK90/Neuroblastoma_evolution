## General settings, such as libraries, samples taken, color codes, etc.
##############################################################################################################################################
## Load libraries that are always needed

library("RColorBrewer")
library(ggplot2); theme_set(theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=6, color="black"),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.75),
                                  axis.ticks = element_line(color = "black", size=0.75),
                                  axis.text = element_text(size=6, color="black")))
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
library(wesanderson)
library(NBevolution)

##############################################################################################################################################
## set directories
data.directory.discovery <- "./Data/discovery/"

snv.directory <- "SNVs/"
indel.directory <- "Indels/"
sv.directory <- "SVs/"
cnv.directory <- "ACEseq/"
custom.script.directory <- "./Analysis_and_plots/"
output.directory <- "./Output/"
meta.data <- "./Meta_data/"
rdata.directory <- "./RData/"
signature.directory <- "./Mutational_signatures/"

if(!dir.exists(rdata.directory)){
  dir.create(rdata.directory)
}

data.directory.validation <- "./Data/validation/"

## directory where the fits for the population genetics model of tumor initiation are stored

fit.directory.initiation <- "./Processed_data/Model_fits_initiation/"

## directory where the fits for the population genetics model of tumor growth are stored

fit.directory.growth <- "./Processed_data/Model_fits_tumor_growth//"


##############################################################################################################################################
## Define samples to use, meta information and colors

sample.information.discovery <- read.xlsx(paste0(meta.data, "Supplementary Table 1.xlsx"), sheet = "Discovery")
rownames(sample.information.discovery) <- sample.information.discovery$Tumor_ID

## classify stage 1, 4S, 2 as low-risk, stage 3 as intermediate risk and stage 4 as high-risk
sample.information.discovery$ManualScore <- "LR"
sample.information.discovery[sample.information.discovery$Stage == "3",]$ManualScore <- "IR"
sample.information.discovery[sample.information.discovery$Stage == "4",]$ManualScore <- "HR"

sample.information.discovery$Stage[sample.information.discovery$Stage %in% c("2.1", "2.2", "2.2000000000000002")] <- "2"

## in rare cases, gains were timed <ECA or between ECA and MRCA; replae <ECA with ECA and ECA<x<MRCA with ECA/MRCA, as we don't distinguish this level for plotting
sample.information.discovery <- replace(sample.information.discovery, sample.information.discovery=="before ECA", "ECA")
sample.information.discovery <- replace(sample.information.discovery, sample.information.discovery=="between ECA and MRCA", "ECA/MRCA")

tumors.discovery <- sample.information.discovery$ID

## Stratify by Lesion type (merge primary and metastasis)
primary.tumors.discovery <- sample.information.discovery$ID[sample.information.discovery$Sample.type %in% c("Primary", "Metastasis")]
relapse.tumors.discovery <- sample.information.discovery$ID[!sample.information.discovery$Sample.type %in% c("Primary", "Metastasis")]

## Stratify by molecular subgroup
telomere.classification.discovery <- sample.information.discovery$Telomere.maintenance.mechanism
names(telomere.classification.discovery) <- rownames(sample.information.discovery)

## Stratify by ploidy
diploid.tumors.discovery <- rownames(sample.information.discovery[sample.information.discovery$Rounded.ploidy==2,])
triploid.tumors.discovery <- rownames(sample.information.discovery[sample.information.discovery$Rounded.ploidy==3,])
tetraploid.tumors.discovery <- rownames(sample.information.discovery[sample.information.discovery$Rounded.ploidy==4,])
diploid.tumors.discovery <- rownames(sample.information.discovery[sample.information.discovery$Rounded.ploidy==2,])


##############################################################################################################################################
## Additional tumors, validation cohort

sample.information.validation <- read.xlsx(paste0(meta.data, "Supplementary Table 1.xlsx"), sheet = "Validation")
rownames(sample.information.validation) <- sample.information.validation$Tumor_ID

tumors.validation <- rownames(sample.information.validation)

primary.tumors.validation <- rownames(sample.information.validation[sample.information.validation$Sample.type %in% c("Primary", "Metastasis"),])

sample.information.validation$ManualScore <- "LR"
sample.information.validation[sample.information.validation$Stage == "3",]$ManualScore <- "IR"
sample.information.validation[sample.information.validation$Stage == "4",]$ManualScore <- "HR"

## in rare cases, gains were timed <ECA or between ECA and MRCA; replae <ECA with ECA and ECA<x<MRCA with ECA/MRCA, as we don't distinguish this level for plotting
sample.information.validation <- replace(sample.information.validation, sample.information.validation=="<ECA", "ECA")
sample.information.validation <- replace(sample.information.validation, sample.information.validation=="ECA<x<MRCA", "ECA/MRCA")

## Stratify by ploidy
diploid.tumors.validation <- rownames(sample.information.validation[sample.information.validation$Rounded.ploidy==2,])
triploid.tumors.validation <- rownames(sample.information.validation[sample.information.validation$Rounded.ploidy==3,])
tetraploid.tumors.validation <- rownames(sample.information.validation[sample.information.validation$Rounded.ploidy==4,])
diploid.tumors.validation <- rownames(sample.information.validation[sample.information.validation$Rounded.ploidy==2,])


## specify the two sample pairs in the data

##  primary tumor of tumor pairs
primary.of.tumor.pairs <- c("NBE11", "NBE51")
## and relapse tumor of tumor pairs
relapse.of.tumor.pairs <- c("NBE66", "NBE78")

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
## Driver genes -- take only genes that have been documented to be NB-relevant

driver.genes <- c("ALK", "ATM", "BRAF", "CCND1", "CDK4", "CDKN2A", "CREBBP", "FGFR1", "LIN28B", "MDM2", "MDM4", "NF1", "PTPN11",
                  "HRAS", "KRAS", "NRAS", "RRAS", "TP53", "MYCN", "TERT", "ATRX")


## genetic positions
gene.positions <- read.delim(paste0(meta.data, "gencode_v19_gene_pos.txt"), header=F, row.names = 1)
driver.genes.with.genetic.pos <- data.frame(GENE=unique(driver.genes), gene.positions[unique(driver.genes),])
colnames(driver.genes.with.genetic.pos) <- c("GENE", "CHROM", "START", "END")

### known non-neutral mutations
known.driver.positions <- read.delim(paste0(meta.data, "Driver_mutations_w_position_DoCM_29_Aug_2020.tsv"),
                                     stringsAsFactors = F, sep="\t")



##############################################################################################################################################
## Based on visual inspection, the purity in tumor NBE40 is not well estimated - adjust it by estimating from the VAF distribution

tumor.to.adjust.purity <- "NBE40"
