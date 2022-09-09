###### Input data to fit neuroblastoma initiation
library(openxlsx)
library(moments)

source("Settings.R")
load(paste0(rdata.directory, "MRCA_timing.RData"))


tmm.tumors <- c(rownames(sample.information.discovery[sample.information.discovery$Telomere.maintenance.mechanism!="None" & 
                                   sample.information.discovery$Sample.type %in% c("Primary", "Metastasis"),]),
  rownames(sample.information.validation[sample.information.validation$Telomere.maintenance.mechanism!="None" & 
                                  sample.information.validation$Sample.type %in% c("Primary", "Metastasis"),]))

no.tmm.tumors <- c(rownames(sample.information.discovery[sample.information.discovery$Telomere.maintenance.mechanism=="None" & 
                                                     sample.information.discovery$Sample.type %in% c("Primary", "Metastasis"),]),
                   rownames(sample.information.validation[sample.information.validation$Telomere.maintenance.mechanism=="None" & 
                                                     sample.information.validation$Sample.type %in% c("Primary", "Metastasis"),]))

## collect the cumulative distribution of MRCA densities
P.MRCA = data.frame(Density=mutation.time.mrca[tmm.tumors,]$Mean)
rownames(P.MRCA) <- tmm.tumors
P.MRCA <- P.MRCA[order(P.MRCA$Density),,drop=F]
P.MRCA <- P.MRCA[!is.na(P.MRCA$Density),,drop=F]
P.MRCA$P <- seq(1, nrow(P.MRCA))/nrow(P.MRCA)
## lower and upper bounds
P.MRCA$P.upper = sapply(P.MRCA$Density, function(x){
  sum(mutation.time.mrca[tmm.tumors,]$Min <= x, na.rm = T)
})/nrow(P.MRCA)
P.MRCA$P.lower = sapply(P.MRCA$Density, function(x){
  sum(mutation.time.mrca[tmm.tumors,]$Max <= x, na.rm = T)
})/nrow(P.MRCA)

## collect the cumulative distribution of ECA densities
P.ECA = data.frame(Density=mutation.time.eca[tmm.tumors,]$Mean)
rownames(P.ECA) <- tmm.tumors
P.ECA <- P.ECA[order(P.ECA$Density),,drop=F]
P.ECA <- P.ECA[!is.na(P.ECA$Density),,drop=F]
## in tetraploid tumors, tetraploidy is likely not the initiating event; thus subset
P.ECA <- P.ECA[setdiff(rownames(P.ECA), c(tetraploid.tumors.discovery, tetraploid.tumors.validation)),,drop=F]
P.ECA$P <- seq(1, nrow(P.ECA))/nrow(P.ECA)
## lower and upper bounds
P.ECA$P.upper = sapply(P.ECA$Density, function(x){
  sum(mutation.time.eca[rownames(P.ECA),]$Min <= x, na.rm = T)
})/nrow(P.ECA)
P.ECA$P.lower = sapply(P.ECA$Density, function(x){
  sum(mutation.time.eca[rownames(P.ECA),]$Max <= x, na.rm = T)
})/nrow(P.ECA)


## collect the cumulative distribution of MRCA densities in low-risk tumors
P.MRCA.lr = data.frame(Density=mutation.time.mrca[no.tmm.tumors,]$Mean)
rownames(P.MRCA.lr) <- no.tmm.tumors
P.MRCA.lr <- P.MRCA.lr[order(P.MRCA.lr$Density),,drop=F]
P.MRCA.lr <- P.MRCA.lr[!is.na(P.MRCA.lr$Density),,drop=F]
P.MRCA.lr$P <- seq(1, nrow(P.MRCA.lr))/nrow(P.MRCA.lr)
## lower and upper bounds
P.MRCA.lr$P.upper = sapply(P.MRCA.lr$Density, function(x){
  sum(mutation.time.mrca[no.tmm.tumors,]$Min <= x, na.rm = T)
})/nrow(P.MRCA.lr)
P.MRCA.lr$P.lower = sapply(P.MRCA.lr$Density, function(x){
  sum(mutation.time.mrca[no.tmm.tumors,]$Max <= x, na.rm = T)
})/nrow(P.MRCA.lr)



save(P.MRCA, P.ECA, P.MRCA.lr, 
     file=paste0(rdata.directory, "Input_data_NB_initiation.RData"))

