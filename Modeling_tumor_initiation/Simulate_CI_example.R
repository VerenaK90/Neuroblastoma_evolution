
## Example simulation

## load fit
fits <- read.csv("./Fitting_results/Expansion_decay.csv")
## and input data
load("./Input_data.RData")
## source model functions
source("./Expansion_decay.R")
## source functions for computing 95% CI
source("./Simulate_CI.R")

parameter.samples <- data.frame(N=fits$par_N, muD1=fits$par_muD1, muD2=fits$par_muD2, mu=fits$par_mu, delta2 = fits$par_delta2,
                                delta1=fits$par_delta1, psurv=fits$par_psurv, r=fits$par_r)
## Simulate incidence curves
sim <- simulate.ci(parameter.samples, NB_origin, NB_origin, NB_origin.eca, mode="expansion_decay")
min.probabilities <- sim[[1]]
max.probabilities <- sim[[2]]
min.probabilities.lr <- sim[[3]]
max.probabilities.lr <- sim[[4]]
min.probabilities.eca <- sim[[5]]
max.probabilities.eca <- sim[[6]]
