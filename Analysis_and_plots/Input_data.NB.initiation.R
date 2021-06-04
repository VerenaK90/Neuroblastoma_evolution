###### Input data to fit neuroblastoma initiation
library(openxlsx)
library(moments)

source(paste0(custom.script.directory, "./Settings.R"))
load(paste0(rdata.directory, "MRCA_timing.RData"))

tmm.tumors <- c(rownames(sample.information.80x[sample.information.80x$Telomere.maintenance.mechanism!="None" & 
                                   sample.information.80x$Location %in% c("Primary", "Metastasis"),]),
  rownames(sample.information.30x[sample.information.30x$Telomere.maintenance.mechanism!="None" & 
                                  sample.information.30x$Location %in% c("Primary", "Metastasis"),]))

no.tmm.tumors <- c(rownames(sample.information.80x[sample.information.80x$Telomere.maintenance.mechanism=="None" & 
                                                     sample.information.80x$Location %in% c("Primary", "Metastasis"),]),
                   rownames(sample.information.30x[sample.information.30x$Telomere.maintenance.mechanism=="None" & 
                                                     sample.information.30x$Location %in% c("Primary", "Metastasis"),]))

NB_origin = sort(mutation.time.mrca[tmm.tumors])
NB_origin.eca = sort(mutation.time.eca[tmm.tumors])
## in tetraploid tumors, tetraploidy is likely not the initiating event
NB_origin.eca = NB_origin.eca[setdiff(names(NB_origin.eca), c(tetraploid.tumors.80x, tetraploid.tumors.30x))]
NB_origin.mrca.lr = sort(mutation.time.mrca[no.tmm.tumors])


###########################################################################

## lower and upper bounds 

NB_origin.lower = sapply(NB_origin, function(x){
  sum(mutation.time.mrca.lower[tmm.tumors] <= x, na.rm = T)
})/length(NB_origin)
NB_origin.upper = sapply(NB_origin, function(x){
  sum(mutation.time.mrca.upper[tmm.tumors] <= x, na.rm = T)
})/length(NB_origin)

NB_origin.eca.lower = sapply(NB_origin.eca, function(x){
  sum(mutation.time.eca.lower[setdiff(names(NB_origin.eca), c(tetraploid.tumors.80x, tetraploid.tumors.30x))] <= x, na.rm = T)
})/length(NB_origin.eca)
NB_origin.eca.upper = sapply(NB_origin.eca, function(x){
  sum(mutation.time.eca.upper[setdiff(names(NB_origin.eca), c(tetraploid.tumors.80x, tetraploid.tumors.30x))] <= x, na.rm = T)
})/length(NB_origin.eca)


NB_origin.mrca.lr.lower = sapply(NB_origin.mrca.lr, function(x){
  sum(mutation.time.mrca.lower[no.tmm.tumors] <= x, na.rm = T)
})/length(NB_origin.mrca.lr)
NB_origin.mrca.lr.upper = sapply(NB_origin.mrca.lr, function(x){
  sum(mutation.time.mrca.upper[no.tmm.tumors] <= x, na.rm = T)
})/length(NB_origin.mrca.lr)


save(NB_origin, NB_origin.lower, NB_origin.upper, 
     NB_origin.eca, NB_origin.eca.lower, NB_origin.eca.upper,
     NB_origin.mrca.lr, NB_origin.mrca.lr.lower, NB_origin.mrca.lr.upper,
     file=paste0(rdata.directory, "Input_data.NB.initiation.RData"))

