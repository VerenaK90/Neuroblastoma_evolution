##Functions to model tumor initiation by 2 events during transient precursor expansion (expansion + contraction)

## Requires the following parameters:
## lambda1: proliferation rate during expansion; set to 1
## delta1: loss rate during expansion
## N: maximal precursor count
## lambda2: proliferation rate during contraction; set to 1
## delta2: loss rate during contraction
## mu1: rate of acquiring 1st oncogenic event
## mu2: rate of acquiring 2nd oncogenic event
## mu: neutral mutation rate
## r: selective advantage of 1st mutant
## p.surv: survival probability of 2nd mutant

##################################################################################
## High-risk tumors: 2 events

## Probability of acquiring both events during the expansion phase (compare Haeno et al., 2007)
## Selective advantage of 1st hit: r (we assume r=1 during expansion as survival probability is anyway high), of 2nd hit: s

P.both.hits.during.expansion <- function(lambda1, delta1, N, t, lambda2, mu1, mu2, p.surv){

  ## stem cell count at time t
  M <- exp((lambda1-delta1)*t)
  ## compute a2, the proliferation rate of the 2nd mutant (lambda1*s) assuming that selection acts on proliferation. Then p.surv = 1-delta1/(s*lambda1) --> s*lambda1 = delta1/(1-p.surv)
  a2 <- delta1/(1-p.surv)
  
  ## require a2 >= 1 (else it would be a selective disadvantage
  if(a2 < 1){
    a2 <- 1
  }
  
  beta <- (1-delta1/lambda1)*mu1/(1-delta1/lambda1)
  
  func <- function(a2, delta1, lambda1){
    integrand <- function(lambda1, delta1, a2, z){
      (1-delta1/a2)/(1-(delta1/a2)*z^((a2-delta1)/(lambda1-delta1)))
    }
    integrate(f = integrand, lower = 0, upper = 1, lambda1=lambda1, delta1=delta1, a2=a2)$value
  }
  
  integrand <- function(M, lambda1, delta1, beta, mu2, x){
    taux <- log(M)/(lambda1-delta1) - log(x)/(lambda1-delta1)
    y <- exp((lambda1-delta1)*taux)
    res <- exp(-beta*(x-1))*(1-exp(-beta))*(1-exp(-mu2*y*func(a2, delta1, lambda1)/(1-delta1/lambda1)))
    res
  }
  if(is.infinite(integrand(M, lambda1, delta1, beta, mu2, 1))){
    return(10^11)
  }
  res <- sum(sapply(seq(1, log10(M)), function(bins){
    integrate(integrand, lower = 10^(bins-1), upper = 10^bins, M=M, lambda1=lambda1, delta1=delta1, beta = beta, mu2=mu2)$value 
  }))
  res
}

## Probability of acquiring the first event during expansion and the second during contraction

P.2nd.driver.1st.hit.during.expansion <- function(lambda1, delta1, N, t, lambda2, delta2, p.surv, mu1, mu2, r){
  
  ## time point at which the population peaks
  t.max <- log(N)/(lambda1-delta1)
  1 - exp(-mu1*mu2*p.surv*lambda1*lambda2/(lambda2-delta2*r)*N*t.max*(exp((lambda2-delta2*r)*t) - 1))
}

## Probability of acquiring both events during contraction

P.2nd.driver.1st.hit.during.decay <- function(lambda1, delta1, N, t, lambda2, delta2, p.surv, mu1, mu2, r){
  
  res <- 1-exp(-mu1*mu2*lambda2^2*p.surv*N/(delta2*(r-1))*(
    1/(lambda2-delta2)*(exp((lambda2-delta2)*t) - 1) - 1/(lambda2 - r*delta2)*(exp((lambda2 - r*delta2)*t)-1)
  ))
  
  return(res)
}


## Joint probability of acquiring both hits at time t

P.2nd.driver <- function(parms, t){
  
  delta1 <- parms$delta1
  N <- 10^parms$N
  delta2 <- parms$delta2
  p.surv <- parms$psurv
  mu1 <- 10^parms$muD1
  mu2 <- 10^parms$muD2
  r <- parms$r
  ## set lambda1, lambda2 to 1
  lambda1 <- 1
  lambda2 <- 1
  
  ## time point at which the population peaks
  t.turning.point <- log(N)/(lambda1-delta1)
  
  ## return Inf if the maximum population size is not reached after 2000 steps
  if(log10(N)/(lambda1-delta1) > 2000){
    return(10^14)
  }
  
  if(t < t.turning.point){
    res <-  P.both.hits.during.expansion(lambda1, delta1, N, t, lambda2, mu1, mu2, p.surv)
  }else{
     res <- 1 -(1-P.both.hits.during.expansion(lambda1, delta1, N, t.turning.point, lambda2, mu1, mu2, p.surv))*
      (1- P.2nd.driver.1st.hit.during.expansion(lambda1, delta1, N, t - t.turning.point, lambda2, delta2, p.surv, mu1, mu2, r))*
      (1-P.2nd.driver.1st.hit.during.decay(lambda1, delta1, N, t - t.turning.point, lambda2, delta2, p.surv, mu1, mu2, r))
  }

  res
}

##################################################################################
## Low-risk tumors: 1 event

P.1st.driver <- function(parms, t){
  
  ## shut off the generation of first-mutants after 2000 mutations (approximately the maximum measured mutation count at MRCA)
  if(t > 2500/parms$mu){
    return(0)
  }
  
  delta1 <- parms$delta1
  N <- 10^parms$N
  delta2 <- parms$delta2
  mu1 <- 10^parms$muD1
  mu2 <- 10^parms$muD2
  r <- parms$r
  lambda1 <- 1
  lambda2 <- 1

  
  ## in case that r*delta2 > lambda2: the first hit cannot expand during involution. Assuming that a sizeable tumor must make up
  ## at least 1% of the tissue, it must've started its growth before N reached 1% of its maximal size. 
  if(lambda2 < delta2*r){
    t.turning.point <- log(N/100)/(lambda1-delta1)
  }else{
    t.turning.point <- log(N)/(lambda1-delta1)
  }
  
  ## return Inf if the population size is not met
  if(log10(N)/(lambda1-delta1) > 2000){
    return(10^14)
  }
  
  ## survival probability of the first mutant during contraction
  p.surv <- 1-delta2*r/lambda2
  if(p.surv < 0){p.surv <- 0}

  if(t < t.turning.point){
    ## survival probability of the first mutant during expansion (same as for normal cells)
    p.surv <- 1-delta1/lambda1
    res <- 1-exp(-mu1*p.surv/(lambda1-delta1)*(exp((lambda1-delta1)*t)-1))
  }else{
    ## survival probability of the first mutant during expansion (same as for normal cells)
    p.surv.1 <- 1-delta1/lambda1
    res <- 1 - exp(-mu1*p.surv.1/(lambda1-delta1)*(exp((lambda1-delta1)*t.turning.point)-1))*
      exp(-mu1*p.surv*N/(lambda2-delta2)*(exp((lambda2-delta2)*(t-t.turning.point)) -1))
  }

  res
  
}

##################################################################################
## Probability of ECA at time t

P.eca.at.t <- function(parms, t1, t2){
  
  delta1 <- parms$delta1
  N <- 10^parms$N
  delta2 <- parms$delta2
  p.surv <- parms$psurv
  mu1 <- 10^parms$muD1
  mu2 <- 10^parms$muD2
  r <- parms$r
  lambda1 <- 1
  lambda2 <- 1
  
  ## time point at which the probability peaks
  t.turning.point <- log(N)/(lambda1-delta1)
  
  ## for cases where t2 > t.turning.point: the number of expected type-1-cells generated during expansion and contraction
  total.count.of.first.hit.during.expansion <- mu1*N*exp((lambda2 - delta2*r)*(t2 - t.turning.point))*t.turning.point
  total.count.of.first.hit.during.decay <-mu1*N/(delta2*(r-1))*(exp((lambda2-delta2)*(t2-t.turning.point)) - exp((lambda2 - delta2*r)*(t2-t.turning.point)))
  
  if(total.count.of.first.hit.during.decay + total.count.of.first.hit.during.expansion ==0){
    return(0)
  }
  if(is.infinite(total.count.of.first.hit.during.decay + total.count.of.first.hit.during.expansion)){
    return(0)
  }
    
  if(t2 < t.turning.point){
    res <- t1/t2
  }else if(t1 <= t.turning.point){
    res <- mu1*N*exp((lambda2-delta2*r)*(t2-t.turning.point))*t1/
      (total.count.of.first.hit.during.expansion + total.count.of.first.hit.during.decay)
  }else if(t1 <= t2){
    res <- mu1*N*exp((lambda2-delta2*r)*(t2-t.turning.point))*t.turning.point/
      (total.count.of.first.hit.during.expansion + total.count.of.first.hit.during.decay) +
      mu1*N/(delta2*(r-1))*(-exp((lambda2 - delta2*r)*(t2-t.turning.point)) + exp(lambda2*(t2-t.turning.point) +
                                                                                  delta2*(t.turning.point + (r-1)*t1 - r*t2)))/
      (total.count.of.first.hit.during.expansion + total.count.of.first.hit.during.decay)
    
  }else{
    res <-0
  }
  res
}


##################################################################################
## Sample the simulated data

## Function to sample the time point of tumor initiation
Sample.timepoint <- function(parms){

  # First, sample a uniformly distributed number; assume total incidence of 10^-5
  x <- runif(1, 0, 1)*10^-5
  ## Then, choose the associated probability
  tmp <- function(x, parms, t){
    P.2nd.driver(parms, t) - x
  }

  if((P.2nd.driver(parms, 10000) >x & P.2nd.driver(parms, 0)  > x)|
     (P.2nd.driver(parms, 10000) <x & P.2nd.driver(parms, 0)  < x)){
      return(100000000)
  }

  t <- uniroot(tmp, c(0, 10000), x=x, parms=parms)$root
  t
}


## Function to sample the number of mutations at tumor initation from a Poisson distribution

Sample.mutations <- function(parms, t){

 n.mutations <- rpois(n=1, lambda=parms$mu*t)
  n.mutations
}


## Function to sample the time point of the ECA
Sample.timepoint.eca <- function(parms, t2){
  # First, sample a uniformly distributed number
  x <- runif(1, 0, 1)

  ## Then, choose the associated probability
  tmp <- function(x, parms, t1, t2){
    P.eca.at.t(parms, t1, t2) - x
  }
  if(P.eca.at.t(parms, t2, t2)==0){
    return(10^8)
  }
  
  t <- uniroot(tmp, c(0, t2), x=x, parms=parms, t2 = t2)$root
  t
}