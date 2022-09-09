##############################################################################################################################################
source("Settings.R")
##############################################################################################################################################
## load data
load(paste0(rdata.directory, "Purity_ploidy.RData"))
load(paste0(rdata.directory, "/Vafs_all_tumors.RData"))
##############################################################################################################################################

tumors.subclonal.evol <- rownames(sample.information.discovery)[sample.information.discovery$Use.to.fit.tumor.expansion=="T"]

for(i in tumors.subclonal.evol){
  print(i)
  
  fits <- list.files(fit.directory.growth, pattern=paste0(i, ".csv"), full.names = T)
  fits <- read.csv(fits)
  
  source(paste0(custom.script.directory, "Neutral_fit_pre_clonal_and_clonal.R"))
  
  if(fit.haploid){
    sim.haploid <- matrix(0, nrow=1000, ncol=length(mySumStatData$haploid))
  }
  if(fit.diploid){
    sim.diploid <- matrix(0, nrow=1000, ncol=length(mySumStatData$diploid))
  }
  if(fit.triploid){
    sim.triploid <- matrix(0, nrow=1000, ncol=length(mySumStatData$triploid))
  }
  if(fit.tetraploid){
    sim.tetraploid <- matrix(0, nrow=1000, ncol=length(mySumStatData$tetraploid))
  }
  
  for(j in 1:nrow(fits)){
    parms <- list(delta=fits$par_delta[j], n_clonal=fits$par_n_clonal[j], mu = fits$par_mu[j])
    output <- myModel(parms)
    
    if(fit.haploid){
      sim.haploid[j,] <- output$haploid
      max.haploid <- max(apply(sim.haploid, 2, max)*haploid.genome.fraction/(3.3*10^9))
    }else{
      max.haploid <- 0
    }
    if(fit.diploid){
      sim.diploid[j,] <- output$diploid
      max.diploid <- max(apply(sim.diploid, 2, max)*diploid.genome.fraction/(3.3*10^9))
    }else{
      max.diploid <- 0
    }
    if(fit.triploid){
      sim.triploid[j,] <- output$triploid
      max.triploid <- max(apply(sim.triploid, 2, max)*triploid.genome.fraction/(3.3*10^9))
    }else{
      max.triploid <- 0
    }
    if(fit.tetraploid){
      sim.tetraploid[j,] <- output$tetraploid
      max.tetraploid <- max(apply(sim.tetraploid, 2, max)*tetraploid.genome.fraction/(3.3*10^9))
    }else{
      max.tetarploid <- 0
    }
    
  }
  
  y.max <- max(max.haploid, max.diploid, max.triploid, max.tetraploid)
  
  p <- list()
  
  if(fit.haploid){
    to.plot <- data.frame(VAF=seq(0.1, 1, 0.05), Data = mySumStatData$haploid*haploid.genome.fraction/(3.3*10^9),
                          Mmin=apply(sim.haploid, 2, quantile, p=0.025)*haploid.genome.fraction/(3.3*10^9),
                          Mmax=apply(sim.haploid, 2, quantile, p=0.975)*haploid.genome.fraction/(3.3*10^9))
    p[[length(p)+1]] <- eval(substitute(ggplot(to.plot, aes(x=VAF, y=Data, 
                                            ymin = Data - sqrt(Data),
                                            ymax = Data + sqrt(Data)
                                            )) + 
      geom_ribbon(aes(ymin=Mmin, 
                      ymax=Mmax),
                  fill="lightslateblue")+ geom_point() + geom_errorbar(width=0.01) + ggtitle(paste("CN=1, weight=", 
                                                                       round(haploid.genome.fraction/(haploid.genome.fraction+diploid.genome.fraction+triploid.genome.fraction+tetraploid.genome.fraction), digits=2))) + 
      scale_y_continuous(limits=c(0, y.max), name="Cumulative number of SNVs"), list(haploid.genome.fraction=haploid.genome.fraction,
                                                                                     diploid.genome.fraction=diploid.genome.fraction,
                                                                                     triploid.genome.fraction=triploid.genome.fraction,
                                                                                     tetraploid.genome.fraction=tetraploid.genome.fraction)))
  }
  if(fit.diploid){
    to.plot <- data.frame(VAF=seq(0.1, 1, 0.05), Data = mySumStatData$diploid*diploid.genome.fraction/(3.3*10^9),
                          Mmin=apply(sim.diploid, 2, quantile, p=0.025)*diploid.genome.fraction/(3.3*10^9),
                          Mmax=apply(sim.diploid, 2, quantile, p=0.975)*diploid.genome.fraction/(3.3*10^9))
    p[[length(p)+1]] <- eval(substitute(ggplot(to.plot, aes(x=VAF, y=Data,   
                                            ymin = Data - sqrt(Data),
                                            ymax = Data + sqrt(Data)
                                            )) + 
      geom_ribbon(aes(ymin=Mmin, 
                      ymax=Mmax),
                  fill="lightslateblue")+ geom_point() + geom_errorbar(width=0.01) + ggtitle(paste("CN=2, weight=", 
                                                                       round(diploid.genome.fraction/(haploid.genome.fraction+diploid.genome.fraction+triploid.genome.fraction+tetraploid.genome.fraction), digits=2))) +
      scale_y_continuous(limits=c(0, y.max), name="Cumulative number of SNVs"),
    list(haploid.genome.fraction=haploid.genome.fraction,
                                                                                   diploid.genome.fraction=diploid.genome.fraction,
                                                                                   triploid.genome.fraction=triploid.genome.fraction,
                                                                                   tetraploid.genome.fraction=tetraploid.genome.fraction)))
  }
  if(fit.triploid){
    to.plot <- data.frame(VAF=seq(0.1, 1, 0.05), Data = mySumStatData$triploid*triploid.genome.fraction/(3.3*10^9),
                          Mmin=apply(sim.triploid, 2, quantile, p=0.025)*triploid.genome.fraction/(3.3*10^9),
                          Mmax=apply(sim.triploid, 2, quantile, p=0.975)*triploid.genome.fraction/(3.3*10^9))
    p[[length(p)+1]] <- eval(substitute(ggplot(to.plot, aes(x=VAF, y=Data, 
                                            ymin = Data - sqrt(Data),
                                            ymax = Data + sqrt(Data)
                                            )) + 
      geom_ribbon(aes(ymin=Mmin, 
                      ymax=Mmax),
                  fill="lightslateblue")+ geom_point() + geom_errorbar(width=0.01)  + ggtitle(paste("CN=3, weight=", 
                                                                        round(triploid.genome.fraction/(haploid.genome.fraction+diploid.genome.fraction+triploid.genome.fraction+tetraploid.genome.fraction), digits=2))) + 
      scale_y_continuous(limits=c(0, y.max), name="Cumulative number of SNVs"), list(haploid.genome.fraction=haploid.genome.fraction,
                                                                                     diploid.genome.fraction=diploid.genome.fraction,
                                                                                     triploid.genome.fraction=triploid.genome.fraction,
                                                                                     tetraploid.genome.fraction=tetraploid.genome.fraction)))
  }
  if(fit.tetraploid){
    to.plot <- data.frame(VAF=seq(0.1, 1, 0.05), Data = mySumStatData$tetraploid*tetraploid.genome.fraction/(3.3*10^9),
                          Mmin=apply(sim.tetraploid, 2, quantile, p=0.025)*tetraploid.genome.fraction/(3.3*10^9),
                          Mmax=apply(sim.tetraploid, 2, quantile, p=0.975)*tetraploid.genome.fraction/(3.3*10^9))
    p[[length(p)+1]] <- eval(substitute(ggplot(to.plot, aes(x=VAF, y=Data, 
                                            ymin = Data - sqrt(Data),
                                            ymax = Data + sqrt(Data)
                                           )) +
      geom_ribbon(aes( ymin=Mmin, 
                       ymax=Mmax),
                  fill="lightslateblue")+ geom_point() + geom_errorbar(width=0.01) + ggtitle(paste("CN=4, weight=", 
                                                                       round(tetraploid.genome.fraction/(haploid.genome.fraction+diploid.genome.fraction+triploid.genome.fraction+tetraploid.genome.fraction), digits=2))) + 
      scale_y_continuous(limits=c(0, y.max), name="Cumulative number of SNVs"),  list(haploid.genome.fraction=haploid.genome.fraction,
                                                                                      diploid.genome.fraction=diploid.genome.fraction,
                                                                                      triploid.genome.fraction=triploid.genome.fraction,
                                                                                      tetraploid.genome.fraction=tetraploid.genome.fraction)))
  }
  
  save(p, file=paste0(data.directory.discovery, i, "/pyABC/Plot_fit.RData"))
}