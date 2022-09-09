#' Simulate mutation densities and their 95% CI using population genetic model of neuroblastoma initiation.
#'
#' This function takes as input the parameter estimates obtained by fitting a population genetics model to NB initiation. It then simulates the model using the parameter estimates, which can be useful for plotting the fitting result.
#' @param parameter.samples a data frame with the parameter samples
#' @param P.MRCA a data frame with the measured densities at MRCA among high-risk tumors (TMM). Rownames must identify the individual tumors. Must contain 1 column named density in units #SNVs/haploid genome
#' @param P.MRCA.lr a data frame with the measured densities at MRCA among low-risk tumors (no TMM). Rownames must identify the individual tumors. Must contain 1 column named density in units #SNVs/haploid genome
#' @param P.ECA  a data frame with the measured densities at ECA among high-risk tumors (TMM). Rownames must identify the individual tumors. Must contain 1 column named density in units #SNVs/haploid genome
#' @param mode Whether tumor initiation should be modeled in a homeostatic tissue ("homeostasis") or in a transient tissue ("decay")
#' @return A list containing lower and upper bounds of the 95% CI of the simulated probabilities at the measured incidences of ECA and MRCA
#' simulate.ci()


simulate.ci <- function(parameter.samples, measured.mutation.times, measured.mutation.times.lr=NA,
                        measured.mutation.times.eca, mode="expansion_only"){

  if(mode=="decay"){
    model <- "contraction"
  }else{
    model <- "homeostasis"
  }
  ## HR tumors MRCA
  simulated.incidence <- matrix(0, nrow=nrow(parameter.samples), ncol=nrow(measured.mutation.times))
  for(i in 1:nrow(parameter.samples)){
    set.seed(Sys.time())
    parms <- as.list(parameter.samples[i,])

    ## rescale measured mutations into time (gamma-distributed)
    time.points <- sapply(sort(measured.mutation.times$Density), function(x){rgamma(n = 1, shape = x, rate = parms$mu)})

    ## then at each time point, take the probability of having acquired the 2 hits
    P <- sapply(time.points, function(t){
      P.2nd.driver(parms, t, model=model)
    })

    simulated.incidence[i,] <- P

  }

  min.probabilities <- apply(simulated.incidence, 2, quantile, p=0.025)
  max.probabilities <- apply(simulated.incidence, 2, quantile, p=0.975)

  ## LR tumors MRCA; only done for the decay scenario
  if(mode=="decay" & !is.na(measured.mutation.times.lr)){
    simulated.incidence <- matrix(0, nrow=nrow(parameter.samples), ncol=length(c(measured.mutation.times.lr$Density,
                                                                                 measured.mutation.times$Density)))
    for(i in 1:nrow(parameter.samples)){
      set.seed(Sys.time())
      parms <- as.list(parameter.samples[i,])

      ## rescale measured mutations into time
      time.points <- sapply(sort(c(measured.mutation.times.lr$Density, measured.mutation.times$Density)), function(x){rgamma(n = 1, shape = x, rate = parms$mu)})
      time.points <- sort(time.points)
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
  simulated.incidence <- matrix(0, nrow=nrow(parameter.samples)*1, ncol=length(measured.mutation.times.eca$Density))
  j <- 1
  for(i in rep(1:nrow(parameter.samples), each=1)){
    set.seed(Sys.time())
    parms <- as.list(parameter.samples[i,])

    ## rescale measured mutations into time; MRCA (gamma distributed)
    time.points <- sapply(sort(measured.mutation.times[rownames(measured.mutation.times.eca),]$Density), function(x){rgamma(n = 1, shape = x, rate = parms$mu)})
    ## and ECA
    time.points.eca <- sapply(sort(measured.mutation.times.eca$Density), function(x){rgamma(n = 1, shape = x, rate = parms$mu)})

    ## then sample a time point for the ECA
    x <- runif(length(time.points), 0, 1)

    ## Then, choose the associated probability
    tmp <- function(x, parms, t1, t2){
      P.eca.at.t(parms, t1, t2, model=model) - x
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

  return(list(min.mrca = min.probabilities,
              max.mrca = max.probabilities,
              min.lr = min.probabilities.lr,
              max.lr = max.probabilities.lr,
              min.eca = min.probabilities.eca,
              max.eca = max.probabilities.eca))

}
