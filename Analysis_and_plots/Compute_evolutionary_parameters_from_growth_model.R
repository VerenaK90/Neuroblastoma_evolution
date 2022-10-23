##############################################################################################################################################
### The fit was done with pyABC. In order to reproduce it, you need to run pyABC using the files
### Neutral_fit_pre_clonal_and_clonal.py, executing Neutral_fit_pre_clonal_and_clonal.R
##############################################################################################################################################

## subset of tumors used to learn the parameters of tumor growth
subset <- sample.information.discovery[sample.information.discovery$Use.to.fit.tumor.expansion=="T" ,]

## from these, extract estimates of the effective mutation rate and subsequently of the loss rate 
deltas <- c()
effective.mutation.rates <- c()
division.rate <- c()
mutation.rate.per.day <- c()
for(i in rownames(subset)){
  
  if(!file.exists(paste0(fit.directory.growth, i, ".csv"))){
    next}
  
  fits <- read.csv(paste0(fit.directory.growth, i, ".csv"))
  
  ## in the growth model, the mutation rate is per 2 cells, but per haploid genome, thus take it as it is (*2/2=1) per cell
  ## The effective mutation rate is per effective division. 
  ## store mean, sd
  effective.mutation.rate <- c(mean(fits$par_mu/(1-fits$par_delta)), sd(fits$par_mu/(1-fits$par_delta)))
  effective.mutation.rates <- cbind(effective.mutation.rates, effective.mutation.rate)
  colnames(effective.mutation.rates)[ncol(effective.mutation.rates)] <- i
  
  ## From this, compute delta; mean and sd
  delta <- 1 - mutation.rate[1]/effective.mutation.rate[1] 
  ## we obtain the error by propagation of uncertainties
  delta[2] <- abs(1/effective.mutation.rate[1] *mutation.rate[2] + mutation.rate[1]/effective.mutation.rate[1]^2 *effective.mutation.rate[2])
  deltas <- cbind(deltas, delta)
  colnames(deltas)[ncol(deltas)] <- i
  
  mutational.burden.at.mrca <- mutation.time.mrca[i,]$Mean
  ## assume that mutation times are roughly normally distributed. Thus the standard deviation would correspond to 1/3.84 of the 95% CI
  mutational.burden.at.mrca[2] <- (mutation.time.mrca[i,]$Max - mutation.time.mrca[i,]$Min)/(2*1.96)
  age <- sample.information.discovery[rownames(sample.information.discovery)==i, "Age"]
  
  n.generations <- mutational.burden.at.mrca*2/mutation.rate[1] + 9*log(10)/(1-delta[1])
  n.generations[2] <- 2/mutation.rate[1]*2*mutational.burden.at.mrca[2] +
    mutational.burden.at.mrca[1]*2/mutation.rate[1]^2*mutation.rate[2] + 9*log(10)/(1-delta[1])^2*delta[2]
  
  ## generations until MRCA 
  n.generations.1 <- mutational.burden.at.mrca[1]*2/mutation.rate[1]
  n.generations.1[2] <- 2/mutation.rate[1]*2*mutational.burden.at.mrca[2]+
    mutational.burden.at.mrca[1]*2/mutation.rate[1]^2*mutation.rate[2] 
  n.generations.2 <- 9*log(10)/(1-delta[1])
  n.generations.2[2] <-  9*log(10)/(1-delta[1])^2*delta[2]
  
  
  ## t.total = age + pregnancy
  t.total <- age + 250
  ## initiation is the number of generations until tumor initiation divided by the total number of generations times the total time
  t.init <- mutational.burden.at.mrca[1]*2/mutation.rate[1]/n.generations[1]*t.total
  t.init[2] <- t.total*(n.generations.1[1]*n.generations.2[1]+n.generations.2[1]*n.generations.1[2])/(n.generations.1[1]+n.generations.2[1])^2
  
  
  division.rate <- cbind(division.rate, c(n.generations[1]/(age + 250), n.generations[2]/(age + 250)))
  colnames(division.rate)[ncol(division.rate)] <- i
  
  ## the mutation rate per day is the product of the division rate and the mutation rate
  mut.per.day <-  division.rate[1,i]*mutation.rate[1]
  mut.per.day[2] <- division.rate[1,i]*mutation.rate[2] + division.rate[2,i]*mutation.rate[1]
  mutation.rate.per.day <-  cbind(mutation.rate.per.day, mut.per.day)
  
}

estimated.mutation.rate.per.day <- c(mean(mutation.rate.per.day[1,subset$Sample.type %in% c("Primary", "Metastasis")]),
                                     sqrt(sum((mutation.rate.per.day[2,subset$Sample.type %in% c("Primary", "Metastasis")])^2))/
                                       sum(subset$Sample.type %in% c("Primary", "Metastasis")))

estimated.mutation.rate.per.day <- c(estimated.mutation.rate.per.day[1] - 1.96*estimated.mutation.rate.per.day[2],
                                     estimated.mutation.rate.per.day[1], 
                                     estimated.mutation.rate.per.day[1] + 1.96*estimated.mutation.rate.per.day[2])


save(estimated.mutation.rate.per.day, file=paste0(rdata.directory, "Estimated_mutation_rate_per_day.RData"))

median(effective.mutation.rates[1,])


## 128
