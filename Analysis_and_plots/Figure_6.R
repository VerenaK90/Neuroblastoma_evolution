## Reproduce Fig. 6
##############################################################################################################################################
### Learn the dynamics of tumor initiation
##############################################################################################################################################
### Load libraries and settings

source("./Settings.R")

load(paste0(rdata.directory, "Estimated_mutation_rate_per_day.RData"))

## source data:
wb <- createWorkbook()
wb.s <- createWorkbook()

## store figure panels
panel.directory <- paste0(output.directory, "Figure6/")

if(!dir.exists(panel.directory)){
  dir.create(panel.directory)
}

##### observed data
## re-run
#source(paste0(custom.script.directory, "Input_data_NB_initiation.R"))
## or load
load(paste0(rdata.directory, "Input_data_NB_initiation.RData"))

mySumStatData <- list(P.MRCA=P.MRCA, P.ECA=P.ECA, P.MRCA.lr=P.MRCA.lr)

load(paste0(rdata.directory, "MRCA_timing.RData"))
load(paste0(rdata.directory,"Clonal_mutations_different_CNs.RData"))
load(paste0(rdata.directory,"Vafs_all_tumors.RData"))

##### take the mutation rate from the fit of tumor initiation
fits <- read.csv(paste0(fit.directory.initiation, "/Expansion_decay_continuous_evol.csv"))
mutation.rate <- c(mean(fits$par_mu*2), sd(fits$par_mu*2))

source(paste0(custom.script.directory, "Compute_evolutionary_parameters_from_growth_model.R"))


##############################################################################################################################################
### The fit was done with pyABC. In order to reproduce it, you need to create the input data by running Input_data.R
### Then run Expansion_decay_continuous_evol.py. You need the files Expansion_decay_continuous_evol.R and Expansion_decay_2_hits_continuous_evol.R
### In analogy for homeostatic fits

##############################################################################################################################################
## Figure 6b Expansion + homeostasis

fits <- read.csv(paste0(fit.directory.initiation, "Expansion_homeostasis_continuous_evol.csv"))

parameter.samples <- data.frame(N=fits$par_N, muD1=fits$par_muD1, muD2=fits$par_muD2, mu=fits$par_mu
                                , delta1=fits$par_delta1, psurv=fits$par_psurv, r=fits$par_r)
## Simulate incidence curves
sim <- simulateCI(parameter.samples = parameter.samples, measured.mutation.times = P.MRCA, measured.mutation.times.eca = P.ECA, mode="homeostasis")
min.probabilities <- sim[[1]]
max.probabilities <- sim[[2]]
min.probabilities.eca <- sim[[5]]
max.probabilities.eca <- sim[[6]]
if(length(which(max.probabilities==1))>1){
  min.probabilities[which(max.probabilities==1)[-1]] <- NA
  max.probabilities[which(max.probabilities==1)[-1]] <- NA
}

to.plot <- data.frame(x = P.MRCA$Density/3.3/10^3,
                      data = P.MRCA$P*10^-5, 
                      sd = (P.MRCA$P.upper - P.MRCA$P.lower)/2/1.95*10^-5,
                      lower = min.probabilities, 
                      upper = max.probabilities,
                      Event = rep("Late", nrow(P.MRCA)))

addWorksheet(wb, "b")
writeData(wb, "b", to.plot)

pdf(paste0(panel.directory, "Figure_6b.pdf"), width=5, height=3)
print(ggplot(to.plot, aes(x=x, y = data, ymin = data-sd, ymax = data +sd, col=Event, fill=Event)) + geom_step() + geom_errorbar()+
        geom_stepribbon(aes(x=x, ymin = lower, ymax = upper), alpha=0.5, col=NA)  +
        scale_fill_manual(values=manual.colors) + scale_color_manual(values=manual.colors) + 
        scale_x_continuous(name="#SSNVs/Mb") + scale_y_continuous(name = "% Neuroblastoma cases") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))+
        coord_cartesian(ylim = c(0, 1.5*10^-5)))

dev.off()

##############################################################################################################################################
## Figure 6d Expansion + decay, ECA

fits <- read.csv(paste0(fit.directory.initiation, "Expansion_decay_continuous_evol.csv"))

parameter.samples <- data.frame(N=fits$par_N, delta1=fits$par_delta1, muD1=fits$par_muD1, 
                                muD2=fits$par_muD2, mu=fits$par_mu,
                                delta2=fits$par_delta2, psurv=fits$par_psurv, r=fits$par_r)


## Simulate incidence curves
sim <- simulateCI(parameter.samples, measured.mutation.times = P.MRCA, measured.mutation.times.eca =  P.ECA, mode="decay")
min.probabilities <- sim[[1]]
max.probabilities <- sim[[2]]

min.probabilities.eca <- sim[[5]]
max.probabilities.eca <- sim[[6]]
if(length(which(max.probabilities==1))>1){
  min.probabilities[which(max.probabilities==1)[-1]] <- NA
  max.probabilities[which(max.probabilities==1)[-1]] <- NA
}

to.plot <- data.frame(x = P.MRCA$Density/3.3/10^3,
                      data = P.MRCA$P*10^-5, 
                      sd = (P.MRCA$P.upper - P.MRCA$P.lower)/2/1.95*10^-5,
                      lower = min.probabilities, 
                      upper = max.probabilities,
                      Event = rep("Late", nrow(P.MRCA)))

to.plot.eca <- data.frame(x = P.ECA$Density/3.3/10^3,
                          data = P.ECA$P, 
                          sd = (P.ECA$P.upper - P.ECA$P.lower)/2/1.95,
                          lower = min.probabilities.eca, 
                          upper = max.probabilities.eca,
                          Event = rep("Early", nrow(P.ECA)))


to.store <- to.plot.eca
to.store$Estimated.weeks <- to.store$x*3.3*10^3*2/7/estimated.mutation.rate.per.day[2] + 2
## get estimate of early hit:
mean(c(mySumStatData$P.ECA$Density)*2/7/estimated.mutation.rate.per.day[2]+2)
sd(c(mySumStatData$P.ECA$Density)*2/7/estimated.mutation.rate.per.day[2]+2)

addWorksheet(wb, "d")
writeData(wb, "d", to.store)


pdf(paste0(panel.directory, "Figure_6d.pdf"), width=5, height=3)

print(ggplot(to.plot.eca, aes(x=x, y = data, ymin = data-sd, ymax = data +sd, col=Event, fill=Event)) + geom_step() + geom_errorbar()+
        geom_stepribbon(aes(x=x, ymin = lower, ymax = upper), alpha=0.5)  +
        scale_fill_manual(values=manual.colors) + scale_color_manual(values=manual.colors) + 
        scale_x_continuous(name="#SSNVs/Mb", limits=c(0, max(to.plot$x)), 
                           sec.axis = sec_axis(~. *3.3*10^3*2/7/estimated.mutation.rate.per.day[2] + 2, name="Estimated weeks p.c.")) + 
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
## Figure 6e Expansion + decay, MRCA


addWorksheet(wb, "e")
to.store <- to.plot
writeData(wb, "e", to.store)

pdf(paste0(panel.directory, "Figure_6e.pdf"), width=5, height=3)
print(ggplot(to.plot, aes(x=x, y = data, ymin = data-sd, ymax = data +sd, col=Event, fill=Event)) + geom_step() + geom_errorbar()+
        geom_stepribbon(aes(x=x, ymin = lower, ymax = upper), alpha=0.5)  +
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
## Figure S5a; parameter estimates

pdf(paste0(panel.directory, "Figure_S5a.pdf"), width=7, height=7, useDingbats = F)

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
  
  if(tmp$xv[1]=="psurv" | tmp$yv[1]=="psurv"){next}
  
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

print(ggarrange(plotlist=p, nrow=10, ncol=10, align="hv"))

dev.off()

##############################################################################################################################################
##### Figure 6f, S5b Plot the retraction of neuroblasts and the exapnsion of the first clone over approximate week of pregnancy

mutation.count <- seq(0, max(P.MRCA$Density))
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

to.store <- N.sim.cells
to.store$Mutations <- to.store$Mutations/3.3/10^3
addWorksheet(wb, "f")
writeData(wb, "f", to.store)


pdf(paste0(panel.directory, "Figure_6f.pdf"), width=7, height=7, useDingbats = F)

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

to.store <- N.sim.cells
to.store$Mutations <- to.store$Mutations/3.3/10^3
addWorksheet(wb.s, "b_Neuroblasts")
writeData(wb.s, "b_Neuroblasts", to.store)

to.store <- M1
addWorksheet(wb.s, "b_M1_cells")
writeData(wb.s, "b_M1_cells", to.store)

pdf(paste0(panel.directory, "Figure_S5b.pdf"), width=7, height=7, useDingbats = F)

p <- ggplot(M1, aes(x=t, ymin=`2.5%`, ymax=`97.5%`)) + geom_ribbon(fill="darkgreen") + 
  geom_ribbon(data=N.sim.cells, aes(x=t, ymin=ymin, ymax=ymax), fill="grey")+
  geom_rect(data=data.frame(xmin = estimated.mutation.rate.per.day[1]*(38+2)*7/2/3.3/10^3, 
                            xmax = estimated.mutation.rate.per.day[3]*(38+2)*7/2/3.3/10^3, ymin = 0, ymax = 1), fill="grey", alpha=0.5,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            inherit.aes = F) + 
  geom_vline(xintercept = estimated.mutation.rate.per.day[2]*(38+2)*7/2/3.3/10^3, col="grey") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous(name="Cell count", sec.axis = sec_axis(~. *3.3*10^3*2/7/estimated.mutation.rate.per.day[2] + 2, name="Estimated weeks p.c.")) 

print(p)

dev.off()




##############################################################################################################################################
## Figure 6g: effective mutation rate


to.plot <- data.frame(Age=subset$Age/365, mu.eff.mean = effective.mutation.rates[1,],
                      Division.rate=division.rate[1,],
                      Subtype=subset$Telomere.maintenance.mechanism,
                      Effective.division.rate=(1-deltas[1,]),
                      Sample.type=subset$Sample.type)

to.plot <- to.plot[!to.plot$Sample.type %in% c("Relapse tumor", "Relapse metastasis"),] 

to.plot <- to.plot[to.plot$Subtype %in% c("MNA", "TERT", "ALT", "None"),]
to.plot$Subtype <- factor(to.plot$Subtype, levels=c("MNA", "TERT", "ALT", "None"))

addWorksheet(wb, "g")
writeData(wb, "g", to.plot) 

pdf(paste0(panel.directory, "Figure_6g.pdf"), width = 5, height=4)

p <- ggplot(to.plot, aes(x=Subtype, y=mu.eff.mean)) + 
  geom_boxplot()+ 
  scale_fill_manual(values=time.colors) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_discrete(name = "") + scale_y_log10()

print(p)


dev.off()

##############################################################################################################################################
## Figure 6h: relative loss rates

pdf(paste0(panel.directory, "Figure_6h.pdf"), width = 5, height=4, useDingbats = F)

## to summarize the data appropriately, I compute for each subclass the standard error of the mean

to.plot <- data.frame(delta.mean = sapply(unique(subset$Telomere.maintenance.mechanism), function(x){
  mean(deltas[1,subset$Telomere.maintenance.mechanism==x & subset$Sample.type %in% c("Primary", "Metastasis")])}), 
  delta.sd =sapply(unique(subset$Telomere.maintenance.mechanism), function(x){
    sqrt(sum((deltas[2,subset$Telomere.maintenance.mechanism==x& subset$Sample.type %in% c("Primary", "Metastasis")])^2))/sum(subset$Telomere.maintenance.mechanism==x& subset$Sample.type %in% c("Primary", "Metastasis"))}),
  Subtype=unique(subset$Telomere.maintenance.mechanism))

to.plot <- to.plot[to.plot$Subtype %in% c("MNA", "TERT", "ALT", "None"),]

to.plot$Subtype <- factor(to.plot$Subtype, levels=c("MNA", "TERT", "ALT", "None"))

addWorksheet(wb, "h")
writeData(wb, "h", to.plot) 


p <- ggplot(to.plot, aes(x=Subtype, y=delta.mean, ymin=delta.mean - delta.sd, ymax=delta.mean+delta.sd)) + 
  geom_pointrange()+ scale_y_continuous(limits=c(0,1)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(p)
dev.off()


## for comparison: tumors without TMM
tumors.no.tmm <- rownames(subset[subset$Telomere.maintenance.mechanism=="None",])
mean(deltas[1,tumors.no.tmm])


##############################################################################################################################################
## Figure 6i: division rate

pdf(paste0(panel.directory, "Figure_6i.pdf"), width = 5, height=4, useDingbats = F)

## to summarize the data appropriately, I compute for each subclass the standard error of the mean


to.plot <- data.frame(division.rate.mean = sapply(unique(subset$Telomere.maintenance.mechanism), function(x){
  mean(division.rate[1,subset$Telomere.maintenance.mechanism==x & subset$Sample.type %in% c("Primary", "Metastasis")])}), 
  division.rate.sd =sapply(unique(subset$Telomere.maintenance.mechanism), function(x){
    sqrt(sum((division.rate[2,subset$Telomere.maintenance.mechanism==x & subset$Sample.type %in% c("Primary", "Metastasis")])^2))/sum(subset$Telomere.maintenance.mechanism==x & subset$Sample.type %in% c("Primary", "Metastasis"))}),
  Telomere.maintenance.mechanism=unique(subset$Telomere.maintenance.mechanism))

to.plot <- to.plot[to.plot$Telomere.maintenance.mechanism %in% c("MNA", "TERT", "ALT", "None"),]

to.plot$Telomere.maintenance.mechanism <- factor(to.plot$Telomere.maintenance.mechanism, levels=c("MNA", "TERT", "ALT", "None"))

addWorksheet(wb, "i")
writeData(wb, "i", to.plot) 

p <- ggplot(to.plot, aes(x=Telomere.maintenance.mechanism, y=division.rate.mean, ymin=division.rate.mean - division.rate.sd,
                         ymax=division.rate.mean+division.rate.sd)) + 
  geom_pointrange()+ scale_y_continuous(limits=c(0, 1.1*max(to.plot$division.rate.mean + to.plot$division.rate.sd)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(p)

dev.off()

##############################################################################################################################################
## Fig. S5c: compare parameter estimates with fits when using clock-like mutations only

fits.clock <- read.csv(paste0(fit.directory.initiation, "Clock_like/Expansion_decay_continuous_evol.csv"))

## compute selective advantage from survival probability
fits.clock$par_s <- fits.clock$par_delta2/(1-fits.clock$par_psurv)
fits.clock$muD1D2 <- sqrt(10^fits.clock$par_muD2*10^fits.clock$par_muD1)
fits.clock$par_muD1 <- 10^fits.clock$par_muD1
fits.clock$par_muD2 <- 10^fits.clock$par_muD2
fits.clock <- fits.clock[,-which(colnames(fits.clock) %in% c("X", "particle_id", "distance", "w"))]

fits.clock. <- as.data.frame(t(hdi(fits.clock, credMass = 0.8)))
fits.clock.$Median <- apply(fits.clock, 2, median)
fits.clock <- fits.clock.
rm(fits.clock.)
fits.clock$Parameter <- rownames(fits.clock)
fits.clock$Mutation <- "Clock"

fits.all <- read.csv(paste0(fit.directory.initiation, "/Expansion_decay_continuous_evol.csv"))

## compute selective advantage from survival probability
fits.all$par_s <- fits.all$par_delta2/(1-fits.all$par_psurv)
fits.all$muD1D2 <- sqrt(10^fits.all$par_muD2*10^fits.all$par_muD1)
fits.all$par_muD1 <- 10^fits.all$par_muD1
fits.all$par_muD2 <- 10^fits.all$par_muD2
fits.all <- fits.all[,-which(colnames(fits.all) %in% c("X", "particle_id", "distance", "w"))]

fits.all. <- as.data.frame(t(hdi(fits.all, credMass = 0.8)))
fits.all.$Median <- apply(fits.all, 2, median)
fits.all <- fits.all.
rm(fits.all.)
fits.all$Parameter <- rownames(fits.all)
fits.all$Mutation <- "All"

fits <- rbind(fits.all, fits.clock)

## parameters to plot
parameters <- c("par_N", "par_delta1", "par_delta2", "par_mu", "par_muD1", "par_muD2", "par_s", "par_r", "muD1D2")

fits <- fits[fits$Parameter %in% parameters,]

p1 <- ggplot(fits[!fits$Parameter %in% c("muD1D2", "par_muD1", "par_muD2"),], 
             aes(x=Mutation, y=Median, ymin=lower, ymax=upper)) + geom_point() + geom_errorbar() +
  facet_wrap(~Parameter, scales="free", ncol=3, nrow=2)

p2 <- ggplot(fits[fits$Parameter %in% c("muD1D2", "par_muD1", "par_muD2"),], 
             aes(x=Mutation, y=Median, ymin=lower, ymax=upper)) + geom_point() + geom_errorbar() + scale_y_continuous(trans = log10_trans()) +
  facet_wrap(~Parameter, scales="free", ncol=3, nrow=1)

pdf(paste0(panel.directory, "Figure_S5c.pdf"), width = 5, height=4, useDingbats = F)

print(ggarrange(p1, p2, nrow=2, heights = c(2,1)))

dev.off()

addWorksheet(wb.s, "c")
writeData(wb.s, "c", fits) 


##############################################################################################################################################
## Fig. S5d: compare real-time estimates of ECA/MRCA between both methods

load(paste0(rdata.directory, "Estimated_mutation_rate_per_day_clock_like.RData"))
estimated.mutation.rate.per.day.clock <- estimated.mutation.rate.per.day
load(paste0(rdata.directory, "Estimated_mutation_rate_per_day.RData"))
estimated.mutation.rate.per.day.all <- estimated.mutation.rate.per.day

load(paste0(rdata.directory, "Input_data_NB_initiation_clock_like.RData"))
P.MRCA.clock <- P.MRCA
P.ECA.clock <- P.ECA
load(paste0(rdata.directory, "Input_data_NB_initiation.RData"))
P.MRCA.all <- P.MRCA
P.ECA.all <- P.ECA

to.plot <- rbind(data.frame(Event = "MRCA",
                      All = P.MRCA.all$Density*2/7/estimated.mutation.rate.per.day.all[2]+2,
                      Clock = P.MRCA.clock$Density*2/7/estimated.mutation.rate.per.day.clock[2]+2),
                 data.frame(Event = "ECA",
                            All = P.ECA.all$Density*2/7/estimated.mutation.rate.per.day.all[2]+2,
                            Clock = P.ECA.clock$Density*2/7/estimated.mutation.rate.per.day.clock[2]+2))

p1 <- ggplot(to.plot[to.plot$Event=="MRCA",], aes(x=All, y=Clock)) + geom_point() +
  scale_x_continuous(name="All SNVs") + scale_y_continuous(name="Clock-like SNVs") +
  ggtitle("MRCA") + geom_abline(slope = 1, intercept = 0, linetype=2) + theme(aspect.ratio = 1)

p2 <- ggplot(to.plot[to.plot$Event=="ECA",], aes(x=All, y=Clock)) + geom_point() +
  scale_x_continuous(name="All SNVs") + scale_y_continuous(name="Clock-like SNVs") +
  ggtitle("ECA") + geom_abline(slope = 1, intercept = 0, linetype=2) + theme(aspect.ratio = 1)

pdf(paste0(panel.directory, "Figure_S5d.pdf"), width = 5, height=4, useDingbats = F)

print(ggarrange(p1, p2, nrow=1))

dev.off()

addWorksheet(wb.s, "d")
writeData(wb.s, "d", to.plot) 


##############################################################################################################################################

saveWorkbook(wb, file = paste0(panel.directory, "Source_data_Fig.6.xlsx"), overwrite=T)

saveWorkbook(wb.s, file = paste0(panel.directory, "Source_data_Fig.S5.xlsx"), overwrite=T)

