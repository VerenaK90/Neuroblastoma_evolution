##### Script to run parameter estimation for NB initiation in a transient population with pyABC
rm(list=ls())

##### load functions

source(paste0(function.directory, "/Model_NB_initiation.R"))

##### observed data

load(rdata.directory, "/Input_data.NB.initiation.RData")

mySumStatData <- list(NB_origin=NB_origin, NB_origin.upper=NB_origin.upper, NB_origin.lower=NB_origin.lower,
                      NB_origin.eca = NB_origin.eca, NB_origin.eca.upper = NB_origin.eca.upper, NB_origin.eca.lower = NB_origin.eca.lower)


##### my model
model <- "homeostasis"

myModel <- function(parms){

   ## sample the same number of tumors as we have data for
   time.points <- sapply(1:length(mySumStatData$NB_origin),function(x){ Sample.timepoint(parms, model=model)})
   ## sample mutation counts at MRCA
   n.mutations.mrca <- sapply(time.points, function(t){Sample.mutations(parms, t)})
   ## sample the time point of the ECA for each MRCA
   time.points.eca <- sapply(time.points[which(names(mySumStatData$NB_origin) %in% names(mySumStatData$NB_origin.eca))], function(t2){
      Sample.timepoint.eca(parms, t2, model=model)
   })
   ## sample the number of mutations at the ECA
   n.mutations.eca <- sapply(time.points.eca, function(t){Sample.mutations(parms, t)})
   
   ## enforce that the incidence overall is approximately 10-5
   incidence.at.age.10 <- P.2nd.driver(parms, 10*365, model=model)
   
   list(n.mutations.mrca=n.mutations.mrca, n.mutations.eca=n.mutations.eca, incidence.at.age.10=incidence.at.age.10)
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
 
    ## Contrast with the measured incidence 
   incidence.mrca <- seq(1/length(simulated.incidence.mrca), 1, length.out=length(simulated.incidence.mrca))
   incidence.eca <- seq(1/length(simulated.incidence.eca), 1, length.out=length(simulated.incidence.eca))

   ## Put more weight on early data points to get the initial shape correctly 
   weights <- replace(rep(1, length(mySumStatData$NB_origin)), mySumStatData$NB_origin/3.3/10^3 <= 0.2, 10)

   ## compute the distance, weighted by the bootstrapped errors
   error.eca <- mySumStatData$NB_origin.eca.upper - mySumStatData$NB_origin.eca.lower
   error.eca[error.eca==0] <- min(error.eca[error.eca!=0])
   distance <- sum(weights*(incidence.mrca - simulated.incidence.mrca)^2/((mySumStatData$NB_origin.upper - mySumStatData$NB_origin.lower)/1.95)^2) +
   sum((incidence.eca - simulated.incidence.eca)^2/((error.eca)/1.95)^2)   + ((modelResult$incidence.at.age.10 - 10^-5)/10^-4)^2

   list(res=distance)
}


## Summary distance: already computed in the distance function
mySummaryDistance <- function(sumStatSample, sumStatData){
  sumStatSample$res

}

