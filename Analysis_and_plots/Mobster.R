##############################################################################################################################################
## load libraries

library(mobster)
library(cowplot)

##############################################################################################################################################
## Input data - discovery set only

#Input data must be a `data.frame` (or, `tibble`) with a columun named `VAF` whose values are numerical $0<x<1$, without `NA` entries. Extra columns will be retained, but not modified.

load(paste0(rdata.directory,"Vafs_all_tumors.RData"))


drivers <- read.delim("NB_drivers_cohort.tsv")

for(i in tumors.discovery){
  
  load(paste0(rdata.directory,"Vafs_all_tumors.RData"))
  load(paste0(rdata.directory,"Clonal_mutations_different_ploidies.RData"))
  
  ploidy.all <- ploidy
  purity.all <- purity
  
  purity <- purity.all[i]
  ploidy <- ploidy.all[i]
  drivers. <- drivers[drivers$SAMPLE==i,,drop=F]
  
  ## merge all mutations in a dataframe. We convert all mutations in haploid mutation counts
  
  dataset <- data.frame(VAF=c())
  
  haploid.genome.fraction <- genome.size.all.tumors[[i]][[1]]
  diploid.genome.fraction <- genome.size.all.tumors[[i]][[2]]
  if(length(genome.size.all.tumors[[i]])>=3){
    triploid.genome.fraction <- genome.size.all.tumors[[i]][[3]]
  }else{
    triploid.genome.fraction <- 0
  }
  if(length(genome.size.all.tumors[[i]])>=4){
    tetraploid.genome.fraction <- genome.size.all.tumors[[i]][[4]]
  }else{
    tetraploid.genome.fraction <- 0
  }
  
  ## require segment length >10^8bp
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
  
  if(length(vafs.all.tumors[[i]][[1]])>1){
    
    vafs.haploid <- vafs.all.tumors[[i]][[1]][,2]/rowSums(vafs.all.tumors[[i]][[1]])
    
    tmp.haploid <- data.frame(GENE=sapply(names(vafs.haploid), function(x){strsplit(x, split="[.]")[[1]][1]}),
                      POS=sapply(names(vafs.haploid), function(x){strsplit(x, split="[.]")[[1]][2]}))
    
    tmp.haploid$is_driver <- apply(tmp.haploid, 1, function(x){
      if(any(drivers.$GENE==x[1] & drivers.$POS==x[2])){
        T
      }else{
        F
      }
    })
    
    tmp.haploid$driver_label <- apply(tmp.haploid, 1, function(x){
      if(any(drivers.$GENE==x[1] & drivers.$POS==x[2])){
        x[1]
      }else{
        "."
      }
    })
    
    vafs.haploid <- vafs.haploid[!is.na(vafs.haploid)]
    depth.haploid <- round(mean(rowSums(vafs.all.tumors[[i]][[1]]), na.rm = T))
    if(length(vafs.haploid)==0){
      fit.haploid <- F
    }else{
      ## convert VAFs into CCFs and then divide by 2 --> pseudo-diploid VAFs
      ## f: cancer cell fraction, rho: purity, then VAF=rho*f/(1*rho +2*(1-rho)), hence f=VAF*(2-rho)/rho and pseudo-diploid VAF=rho*f/2=VAF*(2-rho)/2
      tmp <- vafs.haploid*(2-purity)/2
      dataset <- rbind(dataset, data.frame(GENE=tmp.haploid[names(vafs.haploid),]$GENE, POS=tmp.haploid[names(vafs.haploid),]$POS, VAF=tmp, is_driver=tmp.haploid[names(vafs.haploid),]$is_driver, driver_label=tmp.haploid[names(vafs.haploid),]$driver_label))
    }
  }
  if(fit.diploid && length(vafs.all.tumors[[i]][[2]])>1){
    vafs.diploid <- vafs.all.tumors[[i]][[2]][,2]/rowSums(vafs.all.tumors[[i]][[2]])
    
    tmp.diploid <- data.frame(GENE=sapply(names(vafs.diploid), function(x){strsplit(x, split="[.]")[[1]][1]}),
                              POS=sapply(names(vafs.diploid), function(x){strsplit(x, split="[.]")[[1]][2]}))
    
    tmp.diploid$is_driver <- apply(tmp.diploid, 1, function(x){
      if(any(drivers.$GENE==x[1] & drivers.$POS==x[2])){
        T
      }else{
        F
      }
    })
    
    tmp.diploid$driver_label <- apply(tmp.diploid, 1, function(x){
      if(any(drivers.$GENE==x[1] & drivers.$POS==x[2])){
        x[1]
      }else{
        "."
      }
    })
    vafs.diploid <- vafs.diploid[!is.na(vafs.diploid)]
    
    ## to quantify correctly, we need to cut off higher order peaks and add them to the lower order peaks (otherwise we may underestimate the time prior to
    ## the MRCA as some mutations are "lost" due to amplification)
    depth.diploid <- round(mean(rowSums(vafs.all.tumors[[i]][[2]]), na.rm = T))
    
    homozygous.mutations <- vafs.diploid[vafs.diploid>qbinom(p=0.99, size = depth.diploid, prob = purity/2)/depth.diploid]
    vafs.diploid <- vafs.diploid[vafs.diploid<=qbinom(p=0.99, size = depth.diploid, prob = purity/2)/depth.diploid]
    
    if(length(homozygous.mutations)>0){
      vafs.diploid <- c(vafs.diploid, rep(homozygous.mutations/2, 2))
    }
    dataset <- rbind(dataset, data.frame(GENE=tmp.diploid[names(vafs.diploid),]$GENE, POS=tmp.diploid[names(vafs.diploid),]$POS, VAF=vafs.diploid, is_driver=tmp.diploid[names(vafs.diploid),]$is_driver, driver_label=tmp.diploid[names(vafs.diploid),]$driver_label))
    
  }
  if(fit.triploid && length(vafs.all.tumors[[i]][[3]])>1){
    
    vafs.triploid <- vafs.all.tumors[[i]][[3]][,2]/rowSums(vafs.all.tumors[[i]][[3]])
    
    tmp.triploid <- data.frame(GENE=sapply(names(vafs.triploid), function(x){strsplit(x, split="[.]")[[1]][1]}),
                              POS=sapply(names(vafs.triploid), function(x){strsplit(x, split="[.]")[[1]][2]}))
    
    tmp.triploid$is_driver <- apply(tmp.triploid, 1, function(x){
      if(any(drivers.$GENE==x[1] & drivers.$POS==x[2])){
        T
      }else{
        F
      }
    })
    
    tmp.triploid$driver_label <- apply(tmp.triploid, 1, function(x){
      if(any(drivers.$GENE==x[1] & drivers.$POS==x[2])){
        x[1]
      }else{
        "."
      }
    })
    
    vafs.triploid <- vafs.triploid[!is.na(vafs.triploid)]
    ## cut off higher order peaks from tri- and tetraploid tumors
    depth.triploid <- round(mean(rowSums(vafs.all.tumors[[i]][[3]]), na.rm = T))
    
    vafs.triploid.two.alleles <- vafs.triploid[vafs.triploid>qbinom(p=0.99, size = depth.triploid, prob = purity/(3*purity + 2*(1-purity)))/depth.triploid &
                                                 vafs.triploid<=qbinom(p=0.99, size = depth.triploid, prob = 2*purity/(3*purity + 2*(1-purity)))/depth.triploid]
    vafs.triploid.three.alleles <- vafs.triploid[vafs.triploid>qbinom(p=0.99, size = depth.triploid, prob = 2*purity/(3*purity + 2*(1-purity)))/depth.triploid]
    vafs.triploid <- vafs.triploid[vafs.triploid<=qbinom(p=0.99, size = depth.triploid, prob = purity/(3*purity + 2*(1-purity)))/depth.triploid]
    if(length(vafs.triploid.two.alleles)>0){
      vafs.triploid <- c(vafs.triploid, rep(vafs.triploid.two.alleles*1/2, 2))
    }
    if(length(vafs.triploid.three.alleles)>0){
      vafs.triploid <- c(vafs.triploid, rep(vafs.triploid.three.alleles*1/3, 3))
    }
    ## convert VAFs into CCFs and then divide by 2 --> pseudo-diploid VAFs
    ## f: cancer cell fraction, rho: purity, pi: copy number, then VAF=rho*pi*f/(pi*rho +2*(1-rho)), hence f=VAF*(2 + (pi-2)*rho)/rho and pseudo-diploid VAF=rho*f/2=VAF*(2 + (pi-2)*rho)/2
    tmp <- vafs.triploid*(2+(3-2)*purity)/2
    dataset <- rbind(dataset, data.frame(GENE=tmp.triploid[names(vafs.triploid),]$GENE, POS=tmp.triploid[names(vafs.triploid),]$POS, VAF=tmp, is_driver=tmp.triploid[names(vafs.triploid),]$is_driver, driver_label=tmp.triploid[names(vafs.triploid),]$driver_label))
  }
  if(fit.tetraploid && length(vafs.all.tumors[[i]][[4]])>1){
    
    vafs.tetraploid <- vafs.all.tumors[[i]][[4]][,2]/rowSums(vafs.all.tumors[[i]][[4]])
    
    tmp.tetraploid <- data.frame(GENE=sapply(names(vafs.tetraploid), function(x){strsplit(x, split="[.]")[[1]][1]}),
                               POS=sapply(names(vafs.tetraploid), function(x){strsplit(x, split="[.]")[[1]][2]}))
    
    tmp.tetraploid$is_driver <- apply(tmp.tetraploid, 1, function(x){
      if(any(drivers.$GENE==x[1] & drivers.$POS==x[2])){
        T
      }else{
        F
      }
    })
    
    tmp.tetraploid$driver_label <- apply(tmp.tetraploid, 1, function(x){
      if(any(drivers.$GENE==x[1] & drivers.$POS==x[2])){
        x[1]
      }else{
        "."
      }
    })
    vafs.tetraploid <- vafs.tetraploid[!is.na(vafs.tetraploid)]
    ## cut off higher order peaks from tri- and tetraploid tumors
    depth.tetraploid <- round(mean(rowSums(vafs.all.tumors[[i]][[4]]), na.rm = T))
    
    vafs.tetraploid.two.alleles <- vafs.tetraploid[vafs.tetraploid>qbinom(p=0.99, size = depth.tetraploid, prob = purity/(4*purity + 2*(1-purity)))/depth.tetraploid &
                                                     vafs.tetraploid<=qbinom(p=0.99, size = depth.tetraploid, prob = 2*purity/(4*purity + 2*(1-purity)))/depth.tetraploid]
    vafs.tetraploid.three.alleles <- vafs.tetraploid[vafs.tetraploid>qbinom(p=0.99, size = depth.tetraploid, prob = 2*purity/(4*purity + 2*(1-purity)))/depth.tetraploid &
                                                       vafs.tetraploid<=qbinom(p=0.99, size = depth.tetraploid, prob = 3*purity/(4*purity + 2*(1-purity)))/depth.tetraploid ]
    vafs.tetraploid.four.alleles <- vafs.tetraploid[vafs.tetraploid>qbinom(p=0.99, size = depth.tetraploid, prob = 3*purity/(4*purity + 2*(1-purity)))/depth.tetraploid  ]
    vafs.tetraploid <- vafs.tetraploid[vafs.tetraploid<=qbinom(p=0.99, size = depth.tetraploid, prob = purity/(4*purity + 2*(1-purity)))/depth.tetraploid]
    if(length(vafs.tetraploid.two.alleles)>0){
      vafs.tetraploid <- c(vafs.tetraploid, rep(vafs.tetraploid.two.alleles/2,2))
    }
    if(length(vafs.tetraploid.three.alleles)>0){
      vafs.tetraploid <- c(vafs.tetraploid, rep(vafs.tetraploid.three.alleles*1/3, 3))
    }
    if(length(vafs.tetraploid.four.alleles)>0){
      vafs.tetraploid <- c(vafs.tetraploid, rep(vafs.tetraploid.four.alleles*1/4, 4))
    }
    ## convert VAFs into CCFs and then divide by 2 --> pseudo-diploid VAFs
    ## f: cancer cell fraction, rho: purity, pi: copy number, then VAF=rho*pi*f/(pi*rho +2*(1-rho)), hence f=VAF*(2 + (pi-2)*rho)/rho and pseudo-diploid VAF=rho*f/2=VAF*(2 + (pi-2)*rho)/2
    tmp <- vafs.tetraploid*(2+(4-2)*purity)/2
    dataset <- rbind(dataset, data.frame(GENE=tmp.tetraploid[names(vafs.tetraploid),]$GENE, POS=tmp.tetraploid[names(vafs.tetraploid),]$POS, VAF=tmp, is_driver=tmp.tetraploid[names(vafs.tetraploid),]$is_driver, driver_label=tmp.tetraploid[names(vafs.tetraploid),]$driver_label))
  }
  
  
  dataset <- dataset[dataset$VAF!=1 & !is.na(dataset$VAF) & dataset$VAF >= 0.1,]
  
## run in fast mode first and consider re-running if not satisfied
  fit <- mobster_fit(dataset, auto_setup = "FAST")
  
  ## plot the best model
  pdf(paste0(data.directory.discovery, i, "/Mobster_fits.pdf"))
  
  p <- plot(fit$best, main=i)
  print(p)

  dev.off()
  
  # Assign mutations with at least 80% of probability mass to their maximum responsibility
  clusters = Clusters(
    fit$best, 
    cutoff_assignment = .8
  )
  
  save(fit, clusters, file=paste0(data.directory.discovery, "/", i, "/Mobster_fits.RData"))

}
