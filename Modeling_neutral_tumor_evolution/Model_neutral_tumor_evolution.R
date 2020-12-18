## Clone size distribution according to a linear birth-death process (Bailey, 1964) if starting with a single cell
## lambda: division rate
## delta: loss rate
## n: clone size at t.max
## mu: mutation rate

alpha <- function(lambda, delta, t){
  delta*(exp((lambda - delta)*t) - 1)/(lambda*exp((lambda - delta)*t) - delta)}
beta <- function(lambda, delta, t){
  lambda*(exp((lambda - delta)*t) - 1)/(lambda*exp((lambda - delta)*t) - delta)
}
p.ab <- function(t, lambda, delta, n, t.max){
  (1-alpha(lambda, delta, t.max - t))*(1-beta(lambda, delta, t.max -t))*beta(lambda, delta, t.max - t)^(n-1)
}

## approximate the cumulative distribution of mutations (the number of mutations present in at least n.min and at most n.max cells)
## n.min: minimal clone isze
## n.max: maximal clone size
## mu: mutation rate
## t: start time
## t.max: end time

approx.integrand <- function(t, mu, lambda, delta, n.min, n.max, t.max){
  mu*lambda*exp((lambda - delta)*t)*p.ab(t, lambda, delta, n.max, t.max)/log(beta(lambda, delta, t.max - t)) - 
    mu*lambda*exp((lambda - delta)*t)*p.ab(t, lambda, delta, n.min, t.max)/log(beta(lambda, delta, t.max - t))
}


## model clonal and subclonal mutations together; account for the copy number
## parms: list of parameters, containing 
##        - n_clonal: # clonal mutations
##        - mu: mutation rate
##        - delta: loss rate
## ploidy: tumor ploidy
## depth: average sequencing depth

VAF.model <- function(parms, ploidy, depth){ ## parms: n clonals, delta, mu
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

