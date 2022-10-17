##### Script to run parameter estimation for NB initiation in a transient population with pyABC
rm(list=ls())

source("Settings.R")
##### load functions

library(NBevolution)

##### observed data

load(paste0(rdata.directory, "Input_data_NB_initiation.RData"))

mySumStatData <- list(P.MRCA=P.MRCA, P.ECA=P.ECA, P.MRCA.lr=P.MRCA.lr)


##### my model
model <- "homeostasis"

myModel <- function(parms){

   ## sample the same number of tumors as we have data for
   time.points <- sapply(1:length(mySumStatData$P.MRCA$Density),function(x){ Sample.timepoint(parms, model=model)})
   ## sample mutation counts at MRCA
   n.mutations.mrca <- sapply(time.points, function(t){Sample.mutations(parms, t)})
   ## sample the time point of the ECA for each MRCA
   time.points.eca <- sapply(time.points[which(rownames(mySumStatData$P.MRCA) %in% rownames(mySumStatData$P.ECA))], function(t2){
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
   simulated.incidence.mrca <- sapply(sort(mySumStatData$P.MRCA$Density), function(t){
    sum(modelResult$n.mutations.mrca <= t)
   })/nrow(mySumStatData$P.MRCA)
   ## ... and ECA
    simulated.incidence.eca <- sapply(sort(mySumStatData$P.ECA$Density), function(t){
      sum(modelResult$n.mutations.eca <= t)
    })/nrow(mySumStatData$P.ECA)
 
    ## Contrast with the measured incidence 
   incidence.mrca <- seq(1/length(simulated.incidence.mrca), 1, length.out=length(simulated.incidence.mrca))
   incidence.eca <- seq(1/length(simulated.incidence.eca), 1, length.out=length(simulated.incidence.eca))

   ## Put more weight on early data points to get the initial shape correctly 
   weights <- replace(rep(1, nrow(mySumStatData$P.MRCA)), mySumStatData$P.MRCA$Density/3.3/10^3 <= 0.2, 10)
   
   ## compute the distance, weighted by the bootstrapped errors
   error.eca <- mySumStatData$P.ECA$P.upper - mySumStatData$P.ECA$P.lower
   error.eca[error.eca==0] <- min(error.eca[error.eca!=0])
   distance <- sum(weights*(incidence.mrca - simulated.incidence.mrca)^2/((mySumStatData$P.MRCA$P.upper - mySumStatData$P.MRCA$P.lower)/1.95)^2) +
      sum((incidence.eca - simulated.incidence.eca)^2/((error.eca)/1.95)^2)  + ((modelResult$incidence.at.age.10 - 10^-5)/10^-4)^2

   list(res=distance)
}


## Summary distance: already computed in the distance function
mySummaryDistance <- function(sumStatSample, sumStatData){
  sumStatSample$res

}

