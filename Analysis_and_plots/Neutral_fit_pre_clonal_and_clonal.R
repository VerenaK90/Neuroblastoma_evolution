###################### Infer paramters of a neutral model from tumor samples

rm(list=ls())


############ The fit involves the subclonal tails of the haploid, diploid, triploid and tetraploid fraction

## set tumor
i=
  

load(paste0(rdata.directory, "./Clonal_mutations_different_ploidies.RData"))
load(paste0(rdata.directory, "./Vafs_all_tumors.RData"))

ploidy.all <- ploidy
purity.all <- purity

purity <- purity.all[i]
ploidy <- ploidy.all[i]


## store the genome fraction that is haploid, diploid, triploid and tetraploid and decide whether to fit all or some of them; exclude sex chromosomes

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


if(length(vafs.all.tumors[[i]][[1]])>1){
  vafs.haploid <- vafs.all.tumors[[i]][[1]][,2]/rowSums(vafs.all.tumors[[i]][[1]])
  vafs.haploid <- vafs.haploid[!is.na(vafs.haploid)]
  depth.haploid <- round(mean(rowSums(vafs.all.tumors[[i]][[1]]), na.rm = T))
  if(length(vafs.haploid)==0){
    fit.haploid <- F
  }
}

if(fit.diploid && length(vafs.all.tumors[[i]][[2]])>1){
  vafs.diploid <- vafs.all.tumors[[i]][[2]][,2]/rowSums(vafs.all.tumors[[i]][[2]])
  vafs.diploid <- vafs.diploid[!is.na(vafs.diploid)]
  
  ## to quantify correctly, we need to cut off higher order peaks and add them to the lower order peaks (otherwise we may underestimate the time prior to
  ## the MRCA as some mutations are "lost" due to amplification)
  depth.diploid <- round(mean(rowSums(vafs.all.tumors[[i]][[2]]), na.rm = T))
  
  homozygous.mutations <- vafs.diploid[vafs.diploid>qbinom(p=0.95, size = depth.diploid, prob = purity/2)/depth.diploid]
  vafs.diploid <- vafs.diploid[vafs.diploid<=qbinom(p=0.95, size = depth.diploid, prob = purity/2)/depth.diploid]
  
  if(length(homozygous.mutations)>0){
    vafs.diploid <- c(vafs.diploid, rep(homozygous.mutations/2, 2))
  }
  
  
}
if(fit.triploid && length(vafs.all.tumors[[i]][[3]])>1){
  vafs.triploid <- vafs.all.tumors[[i]][[3]][,2]/rowSums(vafs.all.tumors[[i]][[3]])
  vafs.triploid <- vafs.triploid[!is.na(vafs.triploid)]
  ## cut off higher order peaks from tri- and tetraploid tumors
  depth.triploid <- round(mean(rowSums(vafs.all.tumors[[i]][[3]]), na.rm = T))
  
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
if(fit.tetraploid && length(vafs.all.tumors[[i]][[4]])>1){
  vafs.tetraploid <- vafs.all.tumors[[i]][[4]][,2]/rowSums(vafs.all.tumors[[i]][[4]])
  vafs.tetraploid <- vafs.tetraploid[!is.na(vafs.tetraploid)]
  ## cut off higher order peaks from tri- and tetraploid tumors
  depth.tetraploid <- round(mean(rowSums(vafs.all.tumors[[i]][[4]]), na.rm = T))
  
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


## input data
mySumStatData <- list(haploid = cum.data.haploid, diploid = cum.data.diploid, triploid = cum.data.triploid, tetraploid = cum.data.tetraploid)

## model
myModel <- function(parms){

  ## Simulate the site frequency spectrum accounting for stochastic expansion of individual clones
  ## We simulate the histogram approximately by using integration rather than summation (OK, since steps are small)
  alpha <- function(lambda, delta, t){
    delta*(exp((lambda - delta)*t) - 1)/(lambda*exp((lambda - delta)*t) - delta)}
  beta <- function(lambda, delta, t){
    lambda*(exp((lambda - delta)*t) - 1)/(lambda*exp((lambda - delta)*t) - delta)
  }
  p.ab <- function(t, lambda, delta, n, t.max){
    (1-alpha(lambda, delta, t.max - t))*(1-beta(lambda, delta, t.max -t))*beta(lambda, delta, t.max - t)^(n-1)
  }
  approx.integrand <- function(t, mu, lambda, delta, n.min, n.max, t.max){
    mu*lambda*exp((lambda - delta)*t)*p.ab(t, lambda, delta, n.max, t.max)/log(beta(lambda, delta, t.max - t)) - 
      mu*lambda*exp((lambda - delta)*t)*p.ab(t, lambda, delta, n.min, t.max)/log(beta(lambda, delta, t.max - t))
  }

  ## model clonal and subclonal mutations together; account for the copy number
  rmodel <- function(parms, ploidy, depth){ ## parms: n clonals, delta, mu
    ## sample coverages from Poisson distribution
    coverages <- rpois(n=round(parms$n_clonal*ploidy), lambda=depth)
    ## clonal peak: adjust # mutations for copy number; assume binomial distribution
    vafs.clonal <- rbinom(n=round(parms$n_clonal*ploidy), size = coverages, prob = purity/(ploidy*purity + 2*(1-purity)))/coverages
  
    ## Subclonal mutations
    ## simulate a histogram of bins between 0.05 and 1 with bin size = 0.05; assume a tumor size of 10^9 cells
    vafs.subclonal <- unlist(sapply(seq(0.05, 0.99, 0.05), function(n.min){
      tmp <- integrate(approx.integrand, mu = parms$mu*ploidy, lambda = 1, delta = parms$delta, n.min = n.min*10^9, n.max = (n.min + 0.05)*10^9, t.max = log(10^9)/(1-parms$delta),
                     lower = 0, upper = log(10^9)/(1-parms$delta))$value
      if(tmp>1){
        res <- rep(n.min, round(tmp))
      }else{
        tmp <- sample(x = c(0,1), size = 1, prob = c(1-tmp, tmp))
        res <- rep(n.min, round(tmp))
      }
      unlist(res)
    }))
    if(any(is.na(vafs.subclonal) | is.infinite(vafs.subclonal))){
      return(NA)
    }
    ## Sample coverages
    coverages <- rpois(n=length(vafs.subclonal), lambda=depth)
    ## Assume binomial distribution
    vafs.subclonal <- rbinom(length(vafs.subclonal), size = coverages, prob = vafs.subclonal*purity/(ploidy*purity + 2*(1-purity)))/coverages
  
    return(c(vafs.clonal, vafs.subclonal))
  }


  ## Sample mutations for each copy number separately
  if(fit.haploid){
    model <- rmodel(parms, ploidy = 1, depth = depth.haploid)
    cum.sim.haploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(model >= x)
    })
    if(any(is.na(cum.sim.haploid))){return(Inf)}
  }else{
    cum.sim.haploid <- 0
  }
  if(fit.diploid){
    model <- rmodel(parms, ploidy = 2, depth = depth.diploid)
    cum.sim.diploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(model >= x)
    })
    if(any(is.na(cum.sim.diploid))){return(Inf)}
  }else{
    cum.sim.diploid <- 0
  }
  if(fit.triploid){
    model <- rmodel(parms, ploidy = 3, depth = depth.triploid)
    cum.sim.triploid <- sapply(seq(0.1, 1, 0.05), function(x){
      sum(model >= x)
    })
    if(any(is.na(cum.sim.triploid))){return(Inf)}
  }else{
    cum.sim.triploid <- 0
  }
  if(fit.tetraploid){
    model <- rmodel(parms, ploidy = 4, depth = depth.tetraploid)
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

