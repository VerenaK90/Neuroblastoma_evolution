source("Settings.R")
library(mmsig)

#####################################################################################################################
## Read in the fits by SigProfiler to the discovery cohort

## read in the mutation table
mutations.per.sample <- read.delim(paste0(signature.directory, "discovery/SBS96/Samples.txt"))

mutations.per.sample$Substitution <- sapply(mutations.per.sample$Mutation.Types, function(x){
  x <- paste0(substr(x, 3, 3), substr(x, 5, 5))
})

## and the refit to SBS v3.1
refit <- read.delim(paste0(signature.directory, "discovery/COSMIC_SBS96_Activities_refit.txt"))

## sort by total number of mutations
refit <- refit[order(rowSums(refit[,-1]), decreasing = T),]

refit <- melt(refit, variable.name = "Signature", value.name = "# Mutations")
refit$Samples <- factor(refit$Samples, levels=unique(refit$Samples))
refit$Signature <- factor(refit$Signature, levels=rev(unique(refit$Signature)))

ggplot(refit, aes(x=Samples, y=`# Mutations`, fill=Signature)) + geom_col()

ggplot(refit, aes(x=Samples, y=`# Mutations`, fill=Signature)) + geom_bar(position="fill", stat="identity") 


#####################################################################################################################
### refit the signatures per sample using MMSig

sig_ref <- read.delim(paste0(meta.data, "COSMIC_v3.1_SBS_GRCh37.txt"), sep="\t")

## maximal contribution of each signature

max.contr <- sapply(unique(refit$Signature), function(y){
    max(unlist(sapply(unique(refit$Samples), function(x){
    refit[refit$Samples==x & refit$Signature==y,]$`# Mutations`/sum(refit[refit$Samples==x,]$`# Mutations`)
  })))
})
names(max.contr) <- unique(refit$Signature)

percent.samples.w.sig <- sapply(unique(refit$Signature), function(y){
  sum(unlist(sapply(unique(refit$Samples), function(x){
    refit[refit$Samples==x & refit$Signature==y,]$`# Mutations` > 0/length(unique(refit$Samples))*100
  })))
})
names(percent.samples.w.sig) <- unique(refit$Signature)

## only take sigs that were found in at least 10$ of the samples and that have a minimal contribution of 5% in at least 1 sample

sigs.to.use <- names(max.contr)[max.contr>0.05 & percent.samples.w.sig >= 10]
#####################################################################################################################
### refit the signatures per sample using MMSig

sig_ref <- sig_ref[,c("Type", sigs.to.use)]

sig.colors <- c(wes_palettes$Darjeeling1, wes_palettes$Darjeeling2)[1:(length(sigs.to.use)+1)]
## 1 more color since we later want to merge clock-like sigs into a single one
names(sig.colors) <- c(sigs.to.use, "Clock")

set.seed(1)

input <- mutations.per.sample
rownames(input) <- mutations.per.sample$Mutation.Types
input <- input[,-c(1, ncol(input))]


sig_ref$Trinucleotide <- sapply(sig_ref$Type, function(x){
  x <- paste0(substr(x, 1, 1), substr(x, 3, 3), substr(x, 7, 7))
})

sig_ref <- sig_ref[,c(1, ncol(sig_ref), 2:(ncol(sig_ref)-1))]

sig_out <- mm_fit_signatures(muts.input=input, 
                             sig.input=sig_ref,
                             input.format = "classes",
                             sample.sigt.profs = NULL, 
                             strandbias = F,
                             bootstrap = F,
                             iterations = 1000,
                             refcheck=TRUE,
                             cos_sim_threshold = 0.01,
                             force_include = c("SBS1", "SBS5"),
                             dbg=FALSE) 


write.table(x = sig_out$estimate, file=paste0(signature.directory, "/discovery/MMSig_all_discovery.tsv"),
            sep="\t", quote=F, row.names = T, col.names = T)

save(sig_out, file=paste0(rdata.directory, "/MMSig_all_discovery.RData"))


#####################################################################################################################
### Stratify by clonality variants.

## get purities/ploidies
load(paste0(rdata.directory, "Purity_ploidy.RData"))

input <- data.frame(sample = c(), chr = c(), pos = c(), ref = c(), alt = c(), Clonality=c())

for(i in tumors.discovery){
  
  print(i)
  
  ## Read in ploidy-purity estimates:
  aceseq <- list.files(paste0(data.directory.discovery, i, "/", cnv.directory, "/"), pattern="comb_pro_extra")[1]
  if(is.na(aceseq)){next}
  purity. <- purity[i]
  ploidy. <- ploidy[i]
  
  ## read in the mutation file
  files <- list.files(paste0(data.directory.discovery, i, "/", snv.directory, "/"), pattern="somatic_snvs_conf_8_to_10")[1]
  
  mutations <- read.vcf(paste0(data.directory.discovery, i, "/", snv.directory, "/", files))
  
  mutations$vcf$VAF <- Extract.info.from.vcf(mutations, info="VAF", type="snvs", mutationcaller="DKFZ")
  mutations$vcf$depth <- Extract.info.from.vcf(mutations, info="depth", type="snvs", mutationcaller="DKFZ")
  
  copy.number.info <- read.delim(file=paste0(data.directory.discovery, i, "/", cnv.directory, "/", aceseq), sep="\t", stringsAsFactors = F)
  ## obtain the coverage ratios for the mutations of interest
  cnv.info.per.mutation <- Extract.copy.number.info.per.SSNV(mutations, copy.number.info)
  
  mutations <- as.data.frame(mutations$vcf)
  
  mutations$cov.ratio <- cnv.info.per.mutation$coverage.ratio
  
  ## classify mutations as early clonal, late clonal or subclonal based on VAF
  
  mutations$Clonality <- apply(mutations, 1, function(x){
    
    CN <- round((as.numeric(x["cov.ratio"])*ploidy[i]*purity.+2*(1-purity[i])*(as.numeric(x["cov.ratio"])-1))/(purity[i]))
    
    if(is.na(x["cov.ratio"])|is.na(x["VAF"])|as.numeric(x["cov.ratio"])==0|round(as.numeric(x["cov.ratio"])*ploidy[i])==0){return(NA)}
    
    if(as.numeric(x["VAF"]) < qbinom(p = 0.05, size = as.numeric(x["depth"]),
                                     prob =  1*purity[i]/(purity[i]*CN + 2*(1-purity[i])))/as.numeric(x["depth"])){
      "SC"
      }else{  
      "C"
    }
  })
  
  
  tmp <- mutations
  tmp <- tmp[tmp$CHROM %in% c(1:22),]
  tmp <- data.frame(sample = rep(i, nrow(tmp)), chr = paste0("chr", tmp$CHROM), pos = tmp$POS, 
                    ref = tmp$REF, alt = tmp$ALT, Clonality=tmp$Clonality)
  input <- rbind(input, tmp)
}

input <- input[!is.na(input$Clonality),]

## transform to context matrix

for(clonality in c("C", "SC")){
 
  input. <- as.data.frame(t(mut.to.sigs.input(input[input$Clonality==clonality,], sample.id="sample", chr="chr", pos = "pos", ref = "ref", alt = "alt")))
 ## types must be in the same order in the input and in the reference. Otherwise the tool returns wrong results!
  
   input. <- input.[sig_ref$Type,]
  
    sig_out <- mm_fit_signatures(muts.input=input., 
                                 sig.input=sig_ref,
                                 input.format = "classes",
                                 sample.sigt.profs = NULL, 
                                 strandbias = F,
                                 bootstrap = F,
                                 iterations = 100, # 1000 iterations recommended for stable results
                                 refcheck=TRUE,
                                 cos_sim_threshold = 0.01,
                                 force_include = c("SBS1", "SBS5"),
                                 dbg=FALSE) 

  write.table(x = sig_out$estimate, file=paste0(signature.directory, paste0("/discovery/MMSig_output_discovery_", clonality, ".tsv")),
              sep="\t", quote=F, row.names = T, col.names = T)
}


#####################################################################################################################
### refit also the signatures in the validation data set
#####################################################################################################################
## Read in the fits by SigProfiler to the validation cohort

## read in the mutation table
mutations.per.sample <- read.delim(paste0(signature.directory, "validation/SBS96/Samples.txt"))

mutations.per.sample$Substitution <- sapply(mutations.per.sample$Mutation.Types, function(x){
  x <- paste0(substr(x, 3, 3), substr(x, 5, 5))
})

## and the refit to SBS v3.1
refit <- read.delim(paste0(signature.directory, "validation/Activities/COSMIC_SBS96_Activities_refit.txt"))

## sort by total number of mutations
refit <- refit[order(rowSums(refit[,-1]), decreasing = T),]

refit <- melt(refit, variable.name = "Signature", value.name = "# Mutations")
refit$Samples <- factor(refit$Samples, levels=unique(refit$Samples))
refit$Signature <- factor(refit$Signature, levels=rev(unique(refit$Signature)))

ggplot(refit, aes(x=Samples, y=`# Mutations`, fill=Signature)) + geom_col()

ggplot(refit, aes(x=Samples, y=`# Mutations`, fill=Signature)) + geom_bar(position="fill", stat="identity") 


#####################################################################################################################
### refit the signatures per sample using MMSig

sig_ref <- read.delim(paste0(meta.data, "COSMIC_v3.1_SBS_GRCh37.txt"), sep="\t")

## maximal contribution of each signature

max.contr <- sapply(unique(refit$Signature), function(y){
  max(unlist(sapply(unique(refit$Samples), function(x){
    refit[refit$Samples==x & refit$Signature==y,]$`# Mutations`/sum(refit[refit$Samples==x,]$`# Mutations`)
  })))
})
names(max.contr) <- unique(refit$Signature)

percent.samples.w.sig <- sapply(unique(refit$Signature), function(y){
  sum(unlist(sapply(unique(refit$Samples), function(x){
    refit[refit$Samples==x & refit$Signature==y,]$`# Mutations` > 0/length(unique(refit$Samples))*100
  })))
})
names(percent.samples.w.sig) <- unique(refit$Signature)

## only take sigs that were found in at least 10$ of the samples and that have a minimal contribution of 5% in at least 1 sample

sigs.to.use <- names(max.contr)[max.contr>0.05 & percent.samples.w.sig >= 10]

#####################################################################################################################
### refit the signatures per sample using MMSig

sig_ref <- sig_ref[,c("Type", sigs.to.use)]

## take same colors as before and add for different Sigs between the cohorts
different.sigs <- setdiff(sigs.to.use, names(sig.colors))
## SBS38 and SBS60

sig.colors <- c(sig.colors, 
                c(wes_palettes$Darjeeling1, wes_palettes$Darjeeling2)[(length(sig.colors) + 1): (length(sig.colors) + length(different.sigs))])
## 1 more color since we later want to merge clock-like sigs into a single one
names(sig.colors)[(length(sig.colors) - length(different.sigs) + 1 ):length(sig.colors)] <- different.sigs

# Bootstrapping large datasets with many iterations can significantly increase runtime. 

set.seed(1)

input <- mutations.per.sample
rownames(input) <- mutations.per.sample$Mutation.Types
input <- input[,-c(1, ncol(input))]


sig_ref$Trinucleotide <- sapply(sig_ref$Type, function(x){
  x <- paste0(substr(x, 1, 1), substr(x, 3, 3), substr(x, 7, 7))
})

sig_ref <- sig_ref[,c(1, ncol(sig_ref), 2:(ncol(sig_ref)-1))]

sig_out <- mm_fit_signatures(muts.input=input, 
                             sig.input=sig_ref,
                             input.format = "classes",
                             sample.sigt.profs = NULL, 
                             strandbias = F,
                             bootstrap = F,
                             iterations = 1000, 
                             refcheck=TRUE,
                             cos_sim_threshold = 0.01,
                             force_include = c("SBS1", "SBS5"),
                             dbg=FALSE) 


write.table(x = sig_out$estimate, file=paste0(signature.directory, "/validation/MMSig_all_validation.tsv"),
            sep="\t", quote=F, row.names = T, col.names = T)

save(sig_out, file=paste0(rdata.directory, "/MMSig_all_validation.RData"))



save(sig.colors, file = paste0(rdata.directory, "Sig_colors.RData"))
#####################################################################################################################
## What is the contribution of clock-like mutations per tumor stratified by all, clonal or subclonal variants?

## all mutations
sigs.discovery <- read.delim(paste0(signature.directory, "/discovery/MMSig_all_discovery.tsv"))
## sum up SBS1, SBS5 and SBS40 to a "clock-like" signature
sigs.discovery$Clock <- sigs.discovery$SBS1 + sigs.discovery$SBS5 + sigs.discovery$SBS40
sigs.discovery$Sample <- rownames(sigs.discovery)

## clonal SNVs only
clonal.sigs.discovery <- read.delim(paste0(signature.directory, "/discovery/MMSig_output_discovery_C.tsv"))
clonal.sigs.discovery$Clock <- clonal.sigs.discovery$SBS1 + clonal.sigs.discovery$SBS5 + clonal.sigs.discovery$SBS40
clonal.sigs.discovery$Sample <- rownames(clonal.sigs.discovery)

## subclonal SNVs only
subclonal.sigs.discovery <- read.delim(paste0(signature.directory, "/discovery/MMSig_output_discovery_SC.tsv"))
subclonal.sigs.discovery$Clock <- subclonal.sigs.discovery$SBS1 + subclonal.sigs.discovery$SBS5 + subclonal.sigs.discovery$SBS40
subclonal.sigs.discovery$Sample <- rownames(subclonal.sigs.discovery)

## validation cohort
sigs.validation <- read.delim(paste0(signature.directory, "/validation/MMSig_all_validation.tsv"))
## sum up SBS1, SBS5 and SBS40 to a "clock-like" signature
sigs.validation$Clock <- sigs.validation$SBS1 + sigs.validation$SBS5 + sigs.validation$SBS40
sigs.validation$Sample <- rownames(sigs.validation)

clock.like.all <- c(sigs.discovery$Clock, sigs.validation$Clock)
names(clock.like.all) <- c(sigs.discovery$Sample, sigs.validation$Sample)

clock.like.clonal <- c(clonal.sigs.discovery$Clock, clonal.sigs.discovery$Clock)
names(clock.like.clonal) <- c(clonal.sigs.discovery$Sample, clonal.sigs.discovery$Sample)

clock.like.subclonal <- c(subclonal.sigs.discovery$Clock, subclonal.sigs.discovery$Clock)
names(clock.like.subclonal) <- c(subclonal.sigs.discovery$Sample, subclonal.sigs.discovery$Sample)

save(clock.like.all, clock.like.all, clock.like.subclonal, file=paste0(rdata.directory, "Clock_like_SBS.RData"))



