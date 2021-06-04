## Reproduce Fig. 5
##############################################################################################################################################
### Learn the dynamics of tumor initiation
##############################################################################################################################################
### Load libraries and settings

source(paste0(custom.script.directory, "./Figure_6.R"))

source(paste0(custom.script.directory, "Plotting_function_CI.R"))

load(paste0(rdata.directory, "Estimated_mutation_rate_per_day.RData"))
## source data:
wb <- createWorkbook()
wb.s <- createWorkbook()

## store figure panels
panel.directory <- paste0(output.directory, "Figure5/")

if(!dir.exists(panel.directory)){
  dir.create(panel.directory)
}

##### observed data
## re-run
#source(paste0(custom.script.directory, "Input_data.NB.initiation.R"))
## or load
load(paste0(rdata.directory, "Input_data.NB.initiation.RData"))

mySumStatData <- list(NB_origin=NB_origin, NB_origin.upper=NB_origin.upper, NB_origin.lower=NB_origin.lower,
                      NB_origin.mrca.lr = NB_origin.mrca.lr, NB_origin.mrca.lr.lower = NB_origin.mrca.lr.lower, NB_origin.mrca.lr.upper = NB_origin.mrca.lr.upper,
                      NB_origin.eca = NB_origin.eca, NB_origin.eca.upper = NB_origin.eca.upper, NB_origin.eca.lower = NB_origin.eca.lower)

##############################################################################################################################################
### The fit was done with pyABC. In order to reproduce it, you need to create the input data by running Input_data.NB.initiation
### Then run Expansion_decay_continuous_evol.py. You need the files Expansion_decay_continuous_evol.R and
### Expansion_decay_2_hits_continuous_evol.R
### In analogy for homeostatic fits


##############################################################################################################################################
## Figure 5a Expansion + decay, ECA

fits <- read.csv(paste0(panel.directory, "Expansion_decay_continuous_evol.csv"))

parameter.samples <- data.frame(N=fits$par_N, delta1=fits$par_delta1, muD1=fits$par_muD1, 
                                muD2=fits$par_muD2, mu=fits$par_mu,
                                delta2=fits$par_delta2, psurv=fits$par_psurv, r=fits$par_r)


## Simulate incidence curves
sim <- simulate.ci(parameter.samples, NB_origin, NB_origin.mrca.lr, NB_origin.eca, mode="decay")
min.probabilities <- sim[[1]]
max.probabilities <- sim[[2]]
min.probabilities.lr <- sim[[3]]
max.probabilities.lr <- sim[[4]]
min.probabilities.eca <- sim[[5]]
max.probabilities.eca <- sim[[6]]
if(length(which(max.probabilities==1))>1){
  min.probabilities[which(max.probabilities==1)[-1]] <- NA
  max.probabilities[which(max.probabilities==1)[-1]] <- NA
}

to.plot <- data.frame(x = NB_origin,
                      data = seq(0, 1, length.out = length(NB_origin))*10^-5, 
                      sd = (NB_origin.lower - NB_origin.upper)/2/1.95*10^-5,
                      lower = min.probabilities, 
                      upper = max.probabilities,
                      Event = rep("Late", length(NB_origin)))

to.plot.eca <- data.frame(x = NB_origin.eca,
                          data = seq(0, 1, length.out = length(NB_origin.eca)), 
                          sd = (NB_origin.eca.lower - NB_origin.eca.upper)/2/1.95,
                          lower = min.probabilities.eca, 
                          upper = max.probabilities.eca,
                          Event = rep("Early", length(NB_origin.eca)))


addWorksheet(wb, "a")
writeData(wb, "a", to.plot.eca)


pdf(paste0(panel.directory, "Figure_5a.pdf"), width=5, height=3)

print(ggplot(to.plot.eca, aes(x=x/(3.3*10^9)*10^6, y = data, ymin = data-sd, ymax = data +sd, col=Event, fill=Event)) + geom_step() + geom_errorbar()+
        geom_stepribbon(aes(x=x/(3.3*10^9)*10^6, ymin = lower, ymax = upper), alpha=0.5)  +
        scale_fill_manual(values=manual.colors) + scale_color_manual(values=manual.colors) + 
        scale_x_continuous(name="#SSNVs/Mb", limits=c(0, max(to.plot$x/3.3/10^3)), sec.axis = sec_axis(~. *3.3*10^3*2/7/estimated.mutation.rate.per.day[2] + 2, name="Estimated weeks p.c.")) + 
        scale_y_continuous(name = "% Neuroblastoma cases") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        geom_vline(xintercept = estimated.mutation.rate.per.day[2]*(38-2)*7/2/3.3/10^3, col="grey") +
        geom_rect(data=data.frame(xmin = estimated.mutation.rate.per.day[1]*(38-2)*7/2/3.3/10^3, 
                                  xmax = estimated.mutation.rate.per.day[3]*(38-2)*7/2/3.3/10^3, ymin = 0, ymax = 1), fill="grey", alpha=0.5,
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                  inherit.aes = F))
dev.off()


##############################################################################################################################################
## Figure 5a Expansion + decay, MRCA


addWorksheet(wb, "b")
writeData(wb, "b", to.plot)


pdf(paste0(panel.directory, "Figure_5b.pdf"), width=5, height=3)
print(ggplot(to.plot, aes(x=x/(3.3*10^9)*10^6, y = data, ymin = data-sd, ymax = data +sd, col=Event, fill=Event)) + geom_step() + geom_errorbar()+
        geom_stepribbon(aes(x=x/(3.3*10^9)*10^6, ymin = lower, ymax = upper), alpha=0.5)  +
        coord_cartesian(ylim=c(0, 2*10^-5)) +
        scale_fill_manual(values=manual.colors) + scale_color_manual(values=manual.colors) + 
        scale_x_continuous(name="#SSNVs/Mb", sec.axis = sec_axis(~. *3.3*10^3*2/7/estimated.mutation.rate.per.day[2] + 2, name="Estimated weeks p.c.")) + 
        scale_y_continuous(name = "% Neuroblastoma cases") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
        geom_vline(xintercept = estimated.mutation.rate.per.day[2]*(38-2)*7/2/3.3/10^3, col="grey") +
  geom_rect(data=data.frame(xmin = estimated.mutation.rate.per.day[1]*(38-2)*7/2/3.3/10^3, 
                            xmax = estimated.mutation.rate.per.day[3]*(38-2)*7/2/3.3/10^3, ymin = 0, ymax = 2*10^-5), fill="grey", alpha=0.5,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            inherit.aes = F))

dev.off()



##############################################################################################################################################
## Figure S4b; parameter estimates

pdf(paste0(panel.directory, "Figure_S4b.pdf"), width=7, height=7, useDingbats = F)

## do all correlations
## compute selective advantage from survival probability
fits$par_s <- fits$par_delta2/(1-fits$par_psurv)
parameter.samples$s <- parameter.samples$delta2/(1-parameter.samples$psurv)
## compute geometric mean of mu1 and mu2
parameter.samples$muD1D2 <- log10(sqrt(10^parameter.samples$muD2*10^parameter.samples$muD1))

## parameters to plot
parameters <- c("par_N", "par_delta1", "par_delta2", "par_mu", "par_muD1", "par_muD2", "par_s", "par_r", "muD1D2")

## specify the variables I want to plot
meas_vars <- colnames(parameter.samples)

## a data frame of all combinations of its arguments

controlTable <- data.frame(expand.grid(meas_vars, meas_vars, stringsAsFactors = F))

## rename the columns
colnames(controlTable) <- c("x", "y")

## add the key column
controlTable <- cbind(data.frame(par_key = paste(controlTable[[1]], controlTable[[2]]), stringsAsFactors = F), controlTable)

## create the new data frame
to.plot <- rowrecs_to_blocks(parameter.samples, controlTable)

## re-arrange with facet_grid
splt <- strsplit(to.plot$par_key, split=" ", fixed=TRUE)
to.plot$xv <- vapply(splt, function(si) si[[1]], character(1))
to.plot$yv <- vapply(splt, function(si) si[[2]], character(1))

to.plot$xv <- factor(as.character(to.plot$xv), meas_vars)
to.plot$yv <- factor(as.character(to.plot$yv), meas_vars)

addWorksheet(wb.s, "a")
writeData(wb.s, "a", to.plot)

## arrange manually

to.plot$xaxis <- F
to.plot$yaxis <- F
to.plot$xaxis[to.plot$yv == to.plot$xv[sqrt(length(unique(to.plot$par_key)))]] <- T
to.plot$yaxis[to.plot$xv==to.plot$xv[1]] <- T
to.plot$topm <- F
to.plot$rightm <- F
to.plot$topm[to.plot$yv == to.plot$xv[1]] <- T
to.plot$rightm[to.plot$xv==to.plot$xv[sqrt(length(unique(to.plot$par_key)))]] <- T

p <- list()

## introduce an artificial top row and right column

for(i in 1:(sqrt(length(unique(to.plot$par_key))))){
  p[[length(p)+1]] <- ggplot(data.frame()) + geom_point()+
    theme_bw() + theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"),
                       panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
  
}

for(i in unique(to.plot$par_key)){
  
  
  tmp <- to.plot[to.plot$par_key==i,]
  
  if(tmp$xv=="psurv" | tmp$yv=="psurv"){next}
  
  if(tmp$xv[1]==tmp$yv[1]){
    p[[length(p)+1]] <- ggplot(tmp, aes(x=x)) + 
      geom_histogram() + scale_x_continuous(name=tmp$xv[1]) + scale_y_continuous(name=tmp$yv[1])+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
      theme(legend.position = "none")
    
  }else{
    p[[length(p)+1]] <- ggplot(tmp, aes(x=x, y=y)) + 
      geom_density_2d_filled(col=NA, contour_var = "ndensity", aes( fill = ..level..)) + 
      scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Greens"))(15)) +
      scale_x_continuous(name=tmp$xv[1]) + scale_y_continuous(name=tmp$yv[1])+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
    
  }
  
  ## top-row and right column: adjust margins differently
  if(tmp$rightm[1] & tmp$topm[1]){
    p[[length(p)]] <-  p[[length(p)]] +  theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"))
  }else if(tmp$rightm[1]){
    p[[length(p)]] <-  p[[length(p)]] +  theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"))
  }else if(tmp$topm[1]){
    p[[length(p)]] <-  p[[length(p)]] +  theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"))
  }else{
    p[[length(p)]] <-  p[[length(p)]] +  theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"))
  }

  if(tmp$xaxis[1]==F){
    p[[length(p)]] <-  p[[length(p)]] +  theme(axis.title.x = element_blank(),
                                               axis.text.x = element_blank())
  }
    
  if(tmp$yaxis[1]==F){
    p[[length(p)]] <-  p[[length(p)]] +  theme(axis.title.y = element_blank(),
                                               axis.text.y = element_blank())
  }
  
  if(tmp$rightm[1]){
    p[[length(p)+1]] <- ggplot(data.frame()) + geom_point()+
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
      theme(plot.margin = unit(c(-10, -10, -10, -10), "pt"),
            legend.position = "none")
    
  }
  
}

ggarrange(plotlist=p, nrow=10, ncol=10, align="hv")

dev.off()

##############################################################################################################################################
##### Figure 5c, S4d Plot the retraction of neuroblasts and the exapnsion of the first clone over approximate week of pregnancy

mutation.count <- seq(0, max(NB_origin))
N.sim <- matrix(0, nrow=1000, ncol=length(mutation.count))

hdi.parameters <- hdi(parameter.samples)

for(i in 1:nrow(parameter.samples)){
  delta1 <- parameter.samples[i, "delta1"]
  delta2 <- parameter.samples[i, "delta2"]
  mu <- parameter.samples[i, "mu"]
  N <- 10^parameter.samples[i, "N"]
  
  
  turning.point <- log(N)/(1-delta1)
  N.sim[i,] <- sapply(mutation.count, function(m){
    ## convert to time
    t <- m/mu
    if(t <=turning.point){
      return(exp((1-delta1)*t))
    }else{
      return(N*exp((1-delta2)*(t-turning.point)))
    }
  })
  
}


N.sim <- N.sim[rowSums(N.sim)!=0,]
test <- apply(N.sim, 2, density)

#N.sim.cells <- data.frame(Mutations=rep(mutation.count, nrow(N.sim)), N=as.vector(t(N.sim)), Sim=rep(1:nrow(N.sim), each=ncol(N.sim)))

N.sim.norm <- t(apply(N.sim, 1, function(x){x/max(x)}))
N.sim.cells <- data.frame(Mutations=mutation.count, ymin=apply(N.sim.norm, 2, quantile, p=0.025), ymax=apply(N.sim.norm, 2, quantile, p=0.975))
N.sim.cells$t <- N.sim.cells$Mutations/3.3/10^3

## When does the population peak?

turning.times<- log(10^parameter.samples$N)/(1-parameter.samples$delta1)*parameter.samples$mu*2/7/estimated.mutation.rate.per.day[2] + 2
median(turning.times)
turning.times.upper <- log(10^parameter.samples$N)/(1-parameter.samples$delta1)*parameter.samples$mu*2/7/estimated.mutation.rate.per.day[3] + 2
median(turning.times.upper)
turning.times.lower <- log(10^parameter.samples$N)/(1-parameter.samples$delta1)*parameter.samples$mu*2/7/estimated.mutation.rate.per.day[1] + 2
median(turning.times.lower)

addWorksheet(wb, "c")
writeData(wb, "c", N.sim.cells)


pdf(paste0(panel.directory, "Figure_5c.pdf"), width=7, height=7, useDingbats = F)

ggplot(N.sim.cells, aes(x=t, ymin=ymin, ymax=ymax)) + geom_ribbon(fill="grey") + 
  geom_rect(data=data.frame(xmin = estimated.mutation.rate.per.day[1]*(38+2)*7/2/3.3/10^3, 
                            xmax = estimated.mutation.rate.per.day[3]*(38+2)*7/2/3.3/10^3, ymin = 0, ymax = 1), fill="grey", alpha=0.5,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            inherit.aes = F) + 
  geom_vline(xintercept = estimated.mutation.rate.per.day[2]*(38+2)*7/2/3.3/10^3, col="grey") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(name="Cell count", sec.axis = sec_axis(~. *3.3*10^3*2/7/estimated.mutation.rate.per.day[2] + 2, name="Estimated weeks p.c.")) 


dev.off()

## Simulate evolution of M1 cells

M1.cells <- function(parms, t){
  delta1 <- parms$delta1
  N <- 10^parms$N
  delta2 <- parms$delta2
  mu1 <- 10^parms$muD1
  r <- parms$r
  lambda1 <- 1
  lambda2 <- 1
  
  turning.point <- log(N)/(1-delta1)
  
  if(t <=  turning.point){
    M1 <- mu1*lambda1*t*exp((lambda1-delta1)*t)
  }else if(t>turning.point){
    M1 <- mu1*lambda1*t*N*exp((lambda2 - delta2*r)*(t-turning.point)) +
      mu1*lambda2/(delta2*(r-1))*N*exp((lambda2-delta2*r)*t)/(exp((lambda2-delta2)*turning.point))*(exp(delta2*(r-1)*t)-
                                                                                                      exp(delta2*(r-1)*turning.point))
  }
  M1
}

M1 <- t(apply(parameter.samples, 1, function(parms){
  parms <- as.list(parms)
  cells <- sapply(mutation.count, function(x){
    M1.cells(parms, x/parms$mu)
  })
}))

M1 <- t(apply(M1,1,function(x){x/max(x)}))
M1 <- apply(M1, 2, quantile, p=c(0.025, 0.975))
M1 <- as.data.frame(t(M1))
M1$t <- mutation.count/3.3/10^3

N.sim.cells$t <- N.sim.cells$Mutations/3.3/10^3

addWorksheet(wb.s, "d_Neuroblasts")
writeData(wb.s, "d_Neuroblasts", N.sim.cells)
addWorksheet(wb.s, "d_M1_cells")
writeData(wb.s, "d_M1_cells", M1)

pdf(paste0(panel.directory, "Figure_S4d.pdf"), width=7, height=7, useDingbats = F)

ggplot(M1, aes(x=t, ymin=`2.5%`, ymax=`97.5%`)) + geom_ribbon(fill="darkgreen") + 
  geom_ribbon(data=N.sim.cells, aes(x=t, ymin=ymin, ymax=ymax), fill="grey")+
  geom_rect(data=data.frame(xmin = estimated.mutation.rate.per.day[1]*(38+2)*7/2/3.3/10^3, 
                            xmax = estimated.mutation.rate.per.day[3]*(38+2)*7/2/3.3/10^3, ymin = 0, ymax = 1), fill="grey", alpha=0.5,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            inherit.aes = F) + 
  geom_vline(xintercept = estimated.mutation.rate.per.day[2]*(38+2)*7/2/3.3/10^3, col="grey") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(name="Cell count", sec.axis = sec_axis(~. *3.3*10^3*2/7/estimated.mutation.rate.per.day[2] + 2, name="Estimated weeks p.c.")) 


dev.off()


##############################################################################################################################################
## Figure S4c Expansion + homeostasis

fits <- read.csv(paste0(panel.directory, "Expansion_homeostasis_continuous_evol.csv"))

parameter.samples <- data.frame(N=fits$par_N, muD1=fits$par_muD1, muD2=fits$par_muD2, mu=fits$par_mu
                                , delta1=fits$par_delta1, psurv=fits$par_psurv, r=fits$par_r)
## Simulate incidence curves
sim <- simulate.ci(parameter.samples, NB_origin, NB_origin.mrca.lr, NB_origin.eca, mode="homeostasis")
min.probabilities <- sim[[1]]
max.probabilities <- sim[[2]]
min.probabilities.lr <- sim[[3]]
max.probabilities.lr <- sim[[4]]
min.probabilities.eca <- sim[[5]]
max.probabilities.eca <- sim[[6]]
if(length(which(max.probabilities==1))>1){
  min.probabilities[which(max.probabilities==1)[-1]] <- NA
  max.probabilities[which(max.probabilities==1)[-1]] <- NA
}

to.plot <- data.frame(x = NB_origin,
                      data = seq(0, 1, length.out = length(NB_origin))*10^-5, 
                      sd = (NB_origin.lower - NB_origin.upper)/2/1.95*10^-5,
                      lower = min.probabilities, 
                      upper = max.probabilities,
                      Event = rep("Late", length(NB_origin)))

addWorksheet(wb.s, "c")
writeData(wb.s, "c", to.plot)

pdf(paste0(panel.directory, "Figure_S4c.pdf"), width=5, height=3)
print(ggplot(to.plot, aes(x=x/(3.3*10^9)*10^6, y = data, ymin = data-sd, ymax = data +sd, col=Event, fill=Event)) + geom_step() + geom_errorbar()+
        geom_stepribbon(aes(x=x/(3.3*10^9)*10^6, ymin = lower, ymax = upper), alpha=0.5, col=NA)  +
        scale_fill_manual(values=manual.colors) + scale_color_manual(values=manual.colors) + 
        scale_x_continuous(name="#SSNVs/Mb") + scale_y_continuous(name = "% Neuroblastoma cases") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
        coord_cartesian(ylim = c(0, 1.5*10^-5)))

dev.off()


##############################################################################################################################################

saveWorkbook(wb, file = paste0(panel.directory, "Source_data_Fig.5.xlsx"), overwrite=T)

saveWorkbook(wb.s, file = paste0(panel.directory, "Source_data_Fig.S4.xlsx"), overwrite=T)
