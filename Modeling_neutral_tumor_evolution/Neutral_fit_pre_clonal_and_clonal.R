###################### Infer parameters of a neutral model of tumor evolution from tumor samples

rm(list=ls())

library(bedr)
library(openxlsx)

############ ############ ############ ############ ############ ############ ############ 
############ The fit involves the subclonal tails of the haploid, diploid, triploid and tetraploid fraction

source("./Model_neutral_tumor_evolution.R")

## specify purity and ploidy of the tumor
purity <- 0.8
ploidy <- 3.14

## Required input: SNVs, Indels, ACEseq output. The input files should be organized per tumor
## For each tumor store SNVs and the ACEseq output in subfolders named SNVs, ACEseq

## fill in data directory
data.directory <- ""
## specify tumor
tumor <- ""

############ ############ ############ ############ ############ ############ 

### Assign each mutation to its  copy number state and stor all VAFs on copy number 1,2,3,4. Exclude sex chromosomes

## each entry corresponds to one copy number state
## store the readcounts for refrence and variant reads
readcounts.this.tumor <- list(0)
## and the length of the genome with each copy number state
genome.size.this.tumor <- list(0)

## Read in copy number information and extract ploidy/purity, as before
aceseq <- list.files(paste0(data.directory, tumor, "/ACEseq/"), pattern="comb_pro_extra")[1]
while(is.na(aceseq)){
  aceseq <- list.files(paste0( data.directory, tumor, "/ACEseq/"), pattern="comb_pro_extra")[1]
}

## read in the mutation file
files <- list.files(paste0(data.directory, tumor, "/SNVs/"), pattern="somatic_snvs_conf_8_to_10")[1]
mutations <- read.vcf(paste0( data.directory, tumor, "/SNVs/", files))
mutations$vcf <- mutations$vcf[!mutations$vcf$CHROM %in% c("X", "Y"), ]


copy.number.info <- read.delim(file=paste0(data.directory, tumor, "/ACEseq/", aceseq), sep="\t", stringsAsFactors = F)
copy.number.info <- copy.number.info[!copy.number.info$X.chromosome %in% c("X", "Y"),]
## obtain the coverage ratios for the mutations of interest
coverage.ratios <- apply(mutations$vcf, 1, function(x){
  x <- unlist(x)
  tmp <- copy.number.info[copy.number.info$X.chromosome==x[1],]
  tmp <- tmp[as.numeric(tmp$start) <= as.numeric(x[2]) & as.numeric(tmp$end) >= as.numeric(x[2]),]
  if(nrow(tmp)>0){
    unlist(tmp[28])
  }else{
    NA
  }
})
## and the B-allele frequencies
bafs <- apply(mutations$vcf, 1, function(x){
  x <- unlist(x)
  tmp <- copy.number.info[copy.number.info$X.chromosome==x[1],]
  tmp <- tmp[as.numeric(tmp$start) <= as.numeric(x[2]) & as.numeric(tmp$end) >= as.numeric(x[2]),]
  if(nrow(tmp)>0){
    unlist(tmp[30])
  }else{
    NA
  }
})
## and the genotype
genotype <- apply(mutations$vcf, 1, function(x){
  x <- unlist(x)
  tmp <- copy.number.info[copy.number.info$X.chromosome==x[1],]
  tmp <- tmp[as.numeric(tmp$start) <= as.numeric(x[2]) & as.numeric(tmp$end) >= as.numeric(x[2]),]
  if(nrow(tmp)>0){
    unlist(tmp[36])
  }else{
    NA
  }
})
## and the suggested copy number
tcn <- apply(mutations$vcf, 1, function(x){
  x <- unlist(x)
  tmp <- copy.number.info[copy.number.info$X.chromosome==x[1],]
  tmp <- tmp[as.numeric(tmp$start) <= as.numeric(x[2]) & as.numeric(tmp$end) >= as.numeric(x[2]),]
  if(nrow(tmp)>0){
    as.numeric(unlist(tmp[37]))
  }else{
    NA
  }
})


readcounts <- t(sapply(mutations$vcf$INFO, function(x){
  x <- strsplit(x, split=";")[[1]][6]
  x <- strsplit(x, split="=")[[1]][2]
  x <- as.numeric(strsplit(x, split=",")[[1]])
  x <- unname(x)
  c(sum(x[c(1,2)]), sum(x[c(3,4)]))
}))
rownames(readcounts) <- paste(mutations$vcf$GENE, mutations$vcf$POS, sep=".")


#######################################################################
## Iterate through all copy number states

copy.number.indicator <- c(1:4)

## Plot separately for each copy number
for(k in copy.number.indicator){
  
  expected.coverage.ratio <- (k*purity + (1-purity)*2)/(ploidy*purity+(1-purity)*2)
  
  readcounts. <- readcounts[((coverage.ratios>(expected.coverage.ratio-0.1) & coverage.ratios<(expected.coverage.ratio+0.1) & !is.na(coverage.ratios)) |
                               (tcn ==k & !is.na(tcn))) ,,drop=F]
  
  
  readcounts.this.tumor[[k]] <- readcounts.
  
  
  genome.size <- sum(as.numeric(copy.number.info[(copy.number.info$tcnMeanRaw>(expected.coverage.ratio-0.1) & copy.number.info$tcnMeanRaw<(expected.coverage.ratio+0.1)) |
                                                   (as.numeric(copy.number.info$TCN)==k & !is.na(as.numeric(copy.number.info$TCN))) ,]$end)-
                       as.numeric(copy.number.info[(copy.number.info$tcnMeanRaw>(expected.coverage.ratio-0.1) & copy.number.info$tcnMeanRaw<(expected.coverage.ratio+0.1)) |
                                                     (as.numeric(copy.number.info$TCN)==k & !is.na(as.numeric(copy.number.info$TCN))),]$start))
  
  genome.size.this.tumor[[k]] <- genome.size
  
  
}
  
  
############ ############ ############ ############ ############ ############ 
## Data preparation for fit

## extract the genome fraction that is haploid, diploid, triploid and tetraploid and decide whether to fit all or some of them; exclude sex chromosomes

haploid.genome.fraction <- genome.size.this.tumor[[1]]
diploid.genome.fraction <- genome.size.this.tumor[[2]]
if(length(genome.size.this.tumor)>=3){
  triploid.genome.fraction <- genome.size.this.tumor[[3]]
}else{
  triploid.genome.fraction <- 0
}
if(length(genome.size.this.tumor)>=4){
  tetraploid.genome.fraction <- genome.size.this.tumor[[4]]
}else{
  tetraploid.genome.fraction <- 0
}

## exclude fractions <10^8 bp
if(haploid.genome.fraction < 10^8){
  fit.haploid <- F
}else{
  fit.haploid <- T
}
if(diploid.genome.fraction < 10^8){
  fit.diploid <- F
}else{
  fit.diploid <- T
}
if(triploid.genome.fraction < 10^8){
  fit.triploid <- F
}else{
  fit.triploid <- T
}
if(tetraploid.genome.fraction < 10^8){
  fit.tetraploid <- F
}else{
  fit.tetraploid <- T
}

## prepare data: compute VAFs, merge clonal peaks
if(length(readcounts.this.tumor[[1]])>1){
  vafs.haploid <- readcounts.this.tumor[[1]][,2]/rowSums(readcounts.this.tumor[[1]])
  vafs.haploid <- vafs.haploid[!is.na(vafs.haploid)]
  depth.haploid <- round(mean(rowSums(readcounts.this.tumor[[1]]), na.rm = T))
  if(length(vafs.haploid)==0){
    fit.haploid <- F
  }
}

if(fit.diploid && length(readcounts.this.tumor[[2]])>1){
  vafs.diploid <- readcounts.this.tumor[[2]][,2]/rowSums(readcounts.this.tumor[[2]])
  vafs.diploid <- vafs.diploid[!is.na(vafs.diploid)]
  
  ## to quantify correctly, we need to cut off higher order peaks and add them to the lower order peaks (otherwise we may underestimate the time prior to
  ## the MRCA as some mutations are "lost" due to amplification)
  depth.diploid <- round(mean(rowSums(readcounts.this.tumor[[2]]), na.rm = T))
  
  homozygous.mutations <- vafs.diploid[vafs.diploid>qbinom(p=0.95, size = depth.diploid, prob = purity/2)/depth.diploid]
  vafs.diploid <- vafs.diploid[vafs.diploid<=qbinom(p=0.95, size = depth.diploid, prob = purity/2)/depth.diploid]
  
  if(length(homozygous.mutations)>0){
    vafs.diploid <- c(vafs.diploid, rep(homozygous.mutations/2, 2))
  }
  
  
}
if(fit.triploid && length(readcounts.this.tumor[[3]])>1){
  vafs.triploid <- readcounts.this.tumor[[3]][,2]/rowSums(readcounts.this.tumor[[3]])
  vafs.triploid <- vafs.triploid[!is.na(vafs.triploid)]
  ## cut off higher order peaks from tri- and tetraploid tumors
  depth.triploid <- round(mean(rowSums(readcounts.this.tumor[[3]]), na.rm = T))
  
  vafs.triploid.two.alleles <- vafs.triploid[vafs.triploid>qbinom(p=0.95, size = depth.triploid, prob = purity/(3*purity + 2*(1-purity)))/depth.triploid &
                                               vafs.triploid<=qbinom(p=0.95, size = depth.triploid, prob = 2*purity/(3*purity + 2*(1-purity)))/depth.triploid]
  vafs.triploid.three.alleles <- vafs.triploid[vafs.triploid>qbinom(p=0.95, size = depth.triploid, prob = 2*purity/(3*purity + 2*(1-purity)))/depth.triploid]
  vafs.triploid <- vafs.triploid[vafs.triploid<=qbinom(p=0.95, size = depth.triploid, prob = purity/(3*purity + 2*(1-purity)))/depth.triploid]
  if(length(vafs.triploid.two.alleles)>0){
    vafs.triploid <- c(vafs.triploid, rep(vafs.triploid.two.alleles*1/2, 2))
  }
  if(length(vafs.triploid.three.alleles)>0){
    vafs.triploid <- c(vafs.triploid, rep(vafs.triploid.three.alleles*1/3, 3))
    
  }
  
}
if(fit.tetraploid && length(readcounts.this.tumor[[4]])>1){
  vafs.tetraploid <- readcounts.this.tumor[[4]][,2]/rowSums(readcounts.this.tumor[[4]])
  vafs.tetraploid <- vafs.tetraploid[!is.na(vafs.tetraploid)]
  ## cut off higher order peaks from tri- and tetraploid tumors
  depth.tetraploid <- round(mean(rowSums(readcounts.this.tumor[[4]]), na.rm = T))
  
  vafs.tetraploid.two.alleles <- vafs.tetraploid[vafs.tetraploid>qbinom(p=0.95, size = depth.tetraploid, prob = purity/(4*purity + 2*(1-purity)))/depth.tetraploid &
                                                   vafs.tetraploid<=qbinom(p=0.95, size = depth.tetraploid, prob = 2*purity/(4*purity + 2*(1-purity)))/depth.tetraploid]
  vafs.tetraploid.three.alleles <- vafs.tetraploid[vafs.tetraploid>qbinom(p=0.95, size = depth.tetraploid, prob = 2*purity/(4*purity + 2*(1-purity)))/depth.tetraploid &
                                                     vafs.tetraploid<=qbinom(p=0.95, size = depth.tetraploid, prob = 3*purity/(4*purity + 2*(1-purity)))/depth.tetraploid ]
  vafs.tetraploid.four.alleles <- vafs.tetraploid[vafs.tetraploid>qbinom(p=0.95, size = depth.tetraploid, prob = 3*purity/(4*purity + 2*(1-purity)))/depth.tetraploid  ]
  
  vafs.tetraploid <- vafs.tetraploid[vafs.tetraploid<=qbinom(p=0.95, size = depth.tetraploid, prob = purity/(4*purity + 2*(1-purity)))/depth.tetraploid]
  if(length(vafs.tetraploid.two.alleles)>0){
    vafs.tetraploid <- c(vafs.tetraploid, rep(vafs.tetraploid.two.alleles/2,2))
  }
  if(length(vafs.tetraploid.three.alleles)>0){
    vafs.tetraploid <- c(vafs.tetraploid, rep(vafs.tetraploid.three.alleles*1/3, 3))
  }
  if(length(vafs.tetraploid.four.alleles)>0){
    vafs.tetraploid <- c(vafs.tetraploid, rep(vafs.tetraploid.four.alleles*1/4, 4))
  }
}



### compute the cumulative distribution functions and do upscaling for the entire genome

  ## haploid fraction
  if(fit.haploid){
    cum.data.haploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(vafs.haploid >= x)
    })*3.3*10^9/haploid.genome.fraction
   }else{
    cum.data.haploid <- 0
   }
  if(fit.diploid){
    cum.data.diploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(vafs.diploid >= x)
    })*3.3*10^9/diploid.genome.fraction
   }else{
    cum.data.diploid <- 0
   }
  if(fit.triploid){
    cum.data.triploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(vafs.triploid >= x)
    })*3.3*10^9/triploid.genome.fraction
   }else{
    cum.data.triploid <- 0
   }
  if(fit.tetraploid){
    cum.data.tetraploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(vafs.tetraploid >= x)
    })*3.3*10^9/tetraploid.genome.fraction
   }else{
    cum.data.tetraploid <- 0
   }

############ ############ ############ ############ ############ ############ 

## input data
mySumStatData <- list(haploid = cum.data.haploid, diploid = cum.data.diploid, triploid = cum.data.triploid, tetraploid = cum.data.tetraploid)

## model
myModel <- function(parms){

  ## Simulate the site frequency spectrum accounting for stochastic expansion of individual clones

  ## Sample mutations for each copy number separately
  if(fit.haploid){
    model <- VAF.model(parms, ploidy = 1, depth = depth.haploid)
    cum.sim.haploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(model >= x)
    })
    if(any(is.na(cum.sim.haploid))){return(Inf)}
  }else{
    cum.sim.haploid <- 0
  }
  if(fit.diploid){
    model <- VAF.model(parms, ploidy = 2, depth = depth.diploid)
    cum.sim.diploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(model >= x)
    })
    if(any(is.na(cum.sim.diploid))){return(Inf)}
  }else{
    cum.sim.diploid <- 0
  }
  if(fit.triploid){
    model <- VAF.model(parms, ploidy = 3, depth = depth.triploid)
    cum.sim.triploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(model >= x)
    })
    if(any(is.na(cum.sim.triploid))){return(Inf)}
  }else{
    cum.sim.triploid <- 0
  }
  if(fit.tetraploid){
    model <- VAF.model(parms, ploidy = 4, depth = depth.tetraploid)
    cum.sim.tetraploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(model >= x)
    })
    if(any(is.na(cum.sim.tetraploid))){return(Inf)}
  }else{
    cum.sim.tetraploid <- 0
  }


 list(haploid = cum.sim.haploid, diploid = cum.sim.diploid, triploid = cum.sim.triploid, tetraploid = cum.sim.tetraploid)
}


## Distance to data

mySummaryDistance <- function(sumStatSample, sumStatData){

  
  ## squared distances as error
sqrt( sum(((sumStatSample$haploid - sumStatData$haploid)*haploid.genome.fraction/(haploid.genome.fraction + diploid.genome.fraction + triploid.genome.fraction + tetraploid.genome.fraction))^2) +
    sum(((sumStatSample$diploid - sumStatData$diploid)*diploid.genome.fraction/(haploid.genome.fraction + diploid.genome.fraction + triploid.genome.fraction + tetraploid.genome.fraction))^2) + 
    sum(((sumStatSample$triploid - sumStatData$triploid)*triploid.genome.fraction/(haploid.genome.fraction + diploid.genome.fraction + triploid.genome.fraction + tetraploid.genome.fraction))^2) +
    sum(((sumStatSample$tetraploid - sumStatData$tetraploid)*tetraploid.genome.fraction/(haploid.genome.fraction + diploid.genome.fraction + triploid.genome.fraction + tetraploid.genome.fraction))^2))/
(length(sumStatSample$haploid )+ length(sumStatSample$diploid) + length(sumStatSample$triploid) + length(sumStatSample$tetraploid) ) 
}

