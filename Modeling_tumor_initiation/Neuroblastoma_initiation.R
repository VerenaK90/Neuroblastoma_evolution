##### Script to run parameter estimation for NB initiation in a transient population with pyABC
rm(list=ls())

##### load functions; define model (expansion + decay or expansion + homoestasis); delete other

stem.cell.population <- c("expansion_decay", "expansion_homeostasis")

if(stem.cell.population=="expansion_decay"){
   source("./Expansion_decay.R")
}else if(stem.cell.population=="expansion_homeostasis"){
   source("./Expansion_homeostasis.R")
}

##### observed data

load("./Input_data.RData")

mySumStatData <- list(NB_origin=NB_origin, NB_origin.upper=NB_origin.upper, NB_origin.lower=NB_origin.lower,
                      NB_origin.eca = NB_origin.eca, NB_origin.eca.upper = NB_origin.eca.upper, NB_origin.eca.lower = NB_origin.eca.lower)


##### model to simulate data

myModel <- function(parms){

   ## sample the same number of tumors as we have data for
   time.points <- sapply(1:length(mySumStatData$NB_origin),function(x){ Sample.timepoint(parms)})
   ## sample the mutation counts at MRCA
   n.mutations.mrca <- sapply(time.points, function(t){Sample.mutations(parms, t)})
   ## sample a time point of the ECA for each MRCA
   time.points.eca <- sapply(time.points, function(t2){
      Sample.timepoint.eca(parms, t2)
   })
   ## sample the number of mutations at the ECA
   n.mutations.eca <- sapply(time.points.eca, function(t){Sample.mutations(parms, t)})
   
  list(n.mutations.mrca=n.mutations.mrca, n.mutations.eca=n.mutations.eca)
}



##### Summary statistics

mySummaryStatistics <- function(modelResult){
   ## compute the simulated cumulative incidence at the data points for MRCA...
   simulated.incidence.mrca <- sapply(sort(mySumStatData$NB_origin), function(t){
    sum(modelResult$n.mutations.mrca <= t)
   })/length(mySumStatData$NB_origin)
   ## ... and ECA
   simulated.incidence.eca <- sapply(sort(mySumStatData$NB_origin.eca), function(t){
       sum(modelResult$n.mutations.eca <= t)
   })/length(mySumStatData$NB_origin.eca)
 
   list(simulated.incidence.mrca=simulated.incidence.mrca,
        simulated.incidence.eca=simulated.incidence.eca)
}


## Summary distance: already computed in the distance function
mySummaryDistance <- function(sumStatSample, sumStatData){
   ## Contrast with the measured incidence 
   incidence.mrca <- seq(1/length(sumStatSample$simulated.incidence.mrca), 1, length.out=length(sumStatSample$simulated.incidence.mrca))
   incidence.eca <- seq(1/length(sumStatSample$simulated.incidence.eca), 1, length.out=length(sumStatSample$simulated.incidence.eca))
   
   ## Put more weight on early data points to get the initial shape correctly 
   weights <- replace(rep(1, length(sumStatData$NB_origin)), sumStatData$NB_origin/3.3/10^3 <= 0.2, 10)
   
   ## compute the distance, weighted by the bootstrapped errors
   error.eca <- sumStatData$NB_origin.eca.upper - sumStatData$NB_origin.eca.lower
   error.eca[error.eca==0] <- min(error.eca[error.eca!=0])
   distance <- sum(weights*(incidence.mrca - sumStatSample$simulated.incidence.mrca)^2/((sumStatData$NB_origin.upper - sumStatData$NB_origin.lower)/1.95)^2) +
      sum((incidence.eca - sumStatSample$simulated.incidence.eca)^2/((error.eca)/1.95)^2)
   
   distance

}

