## Compute 95% confidence intervals for simulated data. 
## parameter.samples: a matrix with 1 parameter sample per row; named as defined for the prior probabilities in the python script (N, muD1, muD2, mu, delta1, psurv, r)
## measured.mutation.times: mutation densities at MRCA for high-risk tumors
## measured.mutation.times.eca: mutation densities at ECA for high-risk tumors
## measured.mutation.times.lr: mutation densiites at MRCA for low-risk tumors
## mode: either expansion_homeostasis or expansion_decay, depending on the initiation model

## output: lower and upper bound for the simulated neuroblastoma incidence at the measured densities at MRCA, ECA of high-risk and MRCA of low-risk tumors at the given parameters

simulate.ci <- function(parameter.samples, measured.mutation.times, measured.mutation.times.lr,
                        measured.mutation.times.eca, mode="expansion_decay"){
  
  ## HR tumors MRCA
  simulated.incidence <- matrix(0, nrow=nrow(parameter.samples), ncol=length(measured.mutation.times))
  for(i in 1:nrow(parameter.samples)){
    set.seed(Sys.time())
    parms <- as.list(parameter.samples[i,])
    
    ## rescale measured mutations into time (gamma-distributed)
    time.points <- sapply(sort(measured.mutation.times), function(x){rgamma(n = 1, shape = x, rate = parms$mu)})
    
    ## then at each time point, take the probability of having acquired the 2 hits
    P <- sapply(time.points, function(t){
      P.2nd.driver(parms, t)
    })

    simulated.incidence[i,] <- P
    
  }

  min.probabilities <- apply(simulated.incidence, 2, quantile, p=0.025)
  max.probabilities <- apply(simulated.incidence, 2, quantile, p=0.975)
  
  ## LR tumors MRCA; only done for the decay scenario
  if(mode=="expansion_decay"){
    simulated.incidence <- matrix(0, nrow=nrow(parameter.samples), ncol=length(c(measured.mutation.times.lr, measured.mutation.times)))
    for(i in 1:nrow(parameter.samples)){
      set.seed(Sys.time())
      parms <- as.list(parameter.samples[i,])
      
      ## rescale measured mutations into time
      time.points <- sapply(sort(c(measured.mutation.times.lr, measured.mutation.times)), function(x){rgamma(n = 1, shape = x, rate = parms$mu)})
      
      ## then at each time point, take the probability of having acquired the 1st hit
      
      P <- sapply(time.points, function(t){
        P.1st.driver(parms, t)
      })
      
      simulated.incidence[i,] <- P/max(P)
      
    }
    
    min.probabilities.lr <- apply(simulated.incidence, 2, quantile, p=0.025)
    max.probabilities.lr <- apply(simulated.incidence, 2, quantile, p=0.975)
    
  }else{
    min.probabilities.lr <- NA
    max.probabilities.lr <- NA
  }

  ## ECA
  simulated.incidence <- matrix(0, nrow=nrow(parameter.samples)*1, ncol=length(measured.mutation.times.eca))
  j <- 1
  for(i in rep(1:nrow(parameter.samples), each=1)){
    set.seed(Sys.time())
    parms <- as.list(parameter.samples[i,])
    
    ## rescale measured mutations into time; MRCA (gamma distributed)
    time.points <- sapply(sort(measured.mutation.times[names(measured.mutation.times.eca)]), function(x){rgamma(n = 1, shape = x, rate = parms$mu)})
    ## and ECA
    time.points.eca <- sapply(sort(measured.mutation.times.eca), function(x){rgamma(n = 1, shape = x, rate = parms$mu)})
    
    ## then sample a time point for the ECA
    x <- runif(length(time.points), 0, 1)
    
    ## Then, choose the associated probability
    tmp <- function(x, parms, t1, t2){
      P.eca.at.t(parms, t1, t2) - x
    }

    t <- apply(rbind(time.points, x), 2, function(x){
      uniroot(tmp, c(0, x[1]), x=x[2], parms=parms, t2 = x[1])$root
    })
    
    incidence <- sapply(sort(time.points.eca), function(x){
      sum(t <= x)
    })/length(time.points.eca)

    
    simulated.incidence[j,] <- incidence
    j <- j + 1
  }
  
  min.probabilities.eca <- apply(simulated.incidence, 2, quantile, p=0.025)
  max.probabilities.eca <- apply(simulated.incidence, 2, quantile, p=0.975)
  
  return(list(min = min.probabilities,
              max = max.probabilities,
              min.lr = min.probabilities.lr,
              max.lr = max.probabilities.lr,
              min.eca = min.probabilities.eca,
              max.eca = max.probabilities.eca))
  
}

