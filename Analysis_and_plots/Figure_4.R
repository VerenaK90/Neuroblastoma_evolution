## Reproduce Fig. 4
##############################################################################################################################################
## load settings and libraries
source("./Settings.R")

#source(paste0(custom.script.directory, "Oncoprint.R"))
#source(paste0(custom.script.directory, "Mutation_densities_at_chromosomal_gains.R"))
load(paste0(rdata.directory, "MRCA_timing.RData"))
#load(paste0(rdata.directory, "Estimated_mutation_rate_per_day.RData"))

## source data:
wb <- createWorkbook()
wb.s <- createWorkbook()

## store figure panels
panel.directory <- paste0(output.directory, "Figure4/")

if(!dir.exists(panel.directory)){
  dir.create(panel.directory)
}

source(paste0(custom.script.directory, "Survival_analysis.R"))

## MRCA cutpoint
cutpoint <- 0.05

##########################################################################################################################################
## Figure 4a, b: Survival curves, overall survival for discovery cohort only

## EFS
addWorksheet(wb, "a")
writeData(wb, "a",  data.frame(Time=EFS.fit$time, n.Risk = EFS.fit$n.risk, 
                                 Category=c(rep(names(EFS.fit$strata)[1], EFS.fit$strata[1]),
                                            rep(names(EFS.fit$strata)[2], EFS.fit$strata[2])),
                                 Censored=EFS.fit$n.censor, Event = EFS.fit$n.event, Survival=EFS.fit$surv,
                                 LowerCI=EFS.fit$lower, UpperCI=EFS.fit$upper))

chars <- capture.output(EFS.fit.stats)

writeData(wb, sheet = "a", chars, startCol = 10)

## OS
addWorksheet(wb, "b")
writeData(wb, "b", data.frame(Time=survival.fit$time, n.Risk = survival.fit$n.risk, 
                                          Category=c(rep(names(survival.fit$strata)[1], survival.fit$strata[1]),
                                                     rep(names(survival.fit$strata)[2], survival.fit$strata[2])),
                                          Censored=survival.fit$n.censor, Event = survival.fit$n.event, Survival=survival.fit$surv,
                                          LowerCI=survival.fit$lower, UpperCI=survival.fit$upper))

chars <- capture.output(survival.fit.stats)

writeData(wb, sheet = "b", chars, startCol = 10)


pdf(paste0(panel.directory,"Figure_4a_b.pdf"), useDingbats = F, width=5, height=5)

ggsurvplot(EFS.fit, data = categorized.by.MRCA, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0,10), break.x.by=5, cex=8, linewidth=0.8)

ggsurvplot(survival.fit, data = categorized.by.MRCA, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0,10), break.x.by=5, cex=8, linewidth=0.8) 

dev.off()

##########################################################################################################################################
## Figure 4c, Mutation density at MRCA in the validation cohort

pdf(paste0(panel.directory, "Figure_4c.pdf"), width=3.5, height=3.5, useDingbats = F)

to.plot <- data.frame(MRCA=mutation.time.mrca[primary.tumors.validation,]$Mean)
to.plot$MRCA.time <- ifelse(to.plot$MRCA/3.3/10^3 < 0.05, "Early", "Late")

addWorksheet(wb, "c")
writeData(wb, "c", data.frame(MRCA=to.plot$MRCA/3.3/10^3,
                              MRCA.time=to.plot$MRCA.time))

ggplot(to.plot, aes(x=MRCA/3.3/10^3, fill=MRCA.time)) + geom_histogram(binwidth=0.02) +
  scale_fill_manual(values=c(Early="dodgerblue", Late= "dodgerblue4")) +
  geom_vline(xintercept = cutpoint, linetype=2) + 
  scale_x_continuous(name="#SSNVs/Mb at MRCA") + scale_y_continuous(name="# Tumors") +
  ggtitle("Validation set")

dev.off()


##########################################################################################################################################
## Figure 4d,e: density distributions for MRCA and ECA - validation set

source(paste0(custom.script.directory, "Oncoprint_validation_cohort.R"))

matrix.early.late <- matrix.early.late.validation
## if ECA was not uniquely identified, take the earliest time point
mutation.time.eca[earliest.mutation.time$Sample,]$Mean <- earliest.mutation.time$Mean
mutation.time.eca[earliest.mutation.time$Sample,]$Min <- earliest.mutation.time$Min
mutation.time.eca[earliest.mutation.time$Sample,]$Max <- earliest.mutation.time$Max
####


max.mutation.time.primary <- max(mutation.time.mrca[rownames(sample.information.validation[sample.information.validation$Sample.type %in% c("Primary", "Relapse"),]),]$Mean,
                                 na.rm=T)

sample.information.validation$Telomere.maintenance.mechanism <- factor(sample.information.validation$Telomere.maintenance.mechanism,
                                                                levels=c("MNA", "TERT", "ALT", "Multiple", "None"))

sample.information.validation$ECA.exists <- as.logical(sample.information.validation$ECA.exists)

##### Primary tumors, with and without treatment, metastases
## Plot early and late MRCA separately; further, stratify late cases w.r.t. number of events

plate <- list()
plate.subset.list <- list()

for(ECA.exists in c(T, F)){
  
  subset=sample.information.validation[ sample.information.validation$Sample.type %in% c("Primary", "Metastasis") &
                                   sample.information.validation$ECA.exists==ECA.exists &
                                   mutation.time.mrca[rownames(sample.information.validation),]$Mean/3.3/10^3 >=cutpoint,,drop=F]
  
  if(nrow(subset)==0){next}
  
  
  subset$Telomere.maintenance.mechanism <- factor(subset$Telomere.maintenance.mechanism,
                                                  levels=c("MNA", "TERT", "ALT", "Multiple", "None"))
  
  subset$MRCAtime <- mutation.time.mrca[rownames(subset),]$Mean
  subset <- subset[order(subset$MRCAtime),]
  
  
  to.plot <- cbind(subset, data.frame(MRCA=mutation.time.mrca[rownames(subset),]$Mean/3.3/10^3,
                                      ECA=mutation.time.eca[rownames(subset),]$Mean/3.3/10^3,
                                      MRCA.upper=mutation.time.mrca[rownames(subset),]$Max/3.3/10^3,
                                      MRCA.lower=mutation.time.mrca[rownames(subset),]$Min/3.3/10^3,
                                      ECA.upper=mutation.time.eca[rownames(subset),]$Max/3.3/10^3,
                                      ECA.lower=mutation.time.eca[rownames(subset),]$Min/3.3/10^3))
  
  
  
  if(ECA.exists){
    panel="e_ECA"
  }else{
    panel="e_no_ECA"
  }
  
  addWorksheet(wb, panel)
  writeData(wb, panel, to.plot[,c("ECA", "ECA.lower", "ECA.upper", "MRCA", "MRCA.lower", "MRCA.upper")])
  
  
  
  plate[[length(plate)+1]] <- ggplot(data = to.plot[order(to.plot$MRCA),],
                                     aes(x=MRCA, y=seq(1/length(MRCA),1,length.out = length(MRCA)),
                                         ymin =  sapply(sort(MRCA), function(x){
                                           sum(MRCA.upper <= x)
                                         })/length(MRCA),
                                         ymax= sapply(sort(MRCA), function(x){
                                           sum(MRCA.lower <= x)
                                         })/length(MRCA)
                                     )) +
    stat_ecdf(col=unname(manual.colors["Late"])) +
    geom_stepribbon(fill=unname(manual.colors["Late"]), alpha=0.5, col=NA)+
    scale_x_continuous(name = "Mutations/Mb",
                       limits=c(0, max.mutation.time.primary/3.3/10^3))+
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(name = "Fraction of tumors") +
    ggtitle(paste("# Cases = ", nrow(subset)))
  
  if(ECA.exists){
    plate[[length(plate)]] <- plate[[length(plate)]] + 
      stat_ecdf(data = to.plot[!is.na(to.plot$ECA),][order(to.plot$ECA[!is.na(to.plot$ECA)]),],
                aes(x=ECA, y=seq(1/length(ECA),1,length.out = length(ECA))
                ),
                col=unname(manual.colors["Early"])) +
      geom_stepribbon(data = to.plot[!is.na(to.plot$ECA),][order(to.plot$ECA[!is.na(to.plot$ECA)]),],
                      aes(x=ECA, y=seq(1/length(ECA),1,length.out = length(ECA)),
                          ymin =  sapply(sort(ECA), function(x){
                            sum(ECA.upper <= x)
                          })/length(ECA),
                          ymax= sapply(sort(ECA), function(x){
                            sum(ECA.lower <= x)
                          })/length(ECA)
                      ),
                      fill=unname(manual.colors["Early"]), alpha=0.5, col=NA) 
  }
  
  ## TMM, Clinical subtype & CNVs
  to.plot. <- subset[,c("MRCAtime", "Telomere.maintenance.mechanism", "ManualScore", "Rounded.ploidy")]
  to.plot.$xmin <-c(0, (to.plot.$MRCAtime/3.3/10^3)[-length(to.plot.$MRCAtime)])
  to.plot.$xmax <- to.plot.$MRCAtime/3.3/10^3
  to.plot.$ymin <- 0
  to.plot.$ymax <- 0.23
  
  to.plot.[,c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2")] <- t(matrix.early.late[c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2"),
                                                                                                 rownames(to.plot.)])
  
  ## distinguish cases that are not dateable from cases that are statistically insignificant from ECA
  for(chr.change in c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2")){
    tmp <- not.dateable[not.dateable$Gain==chr.change,]
    to.plot.[intersect(rownames(to.plot.), tmp$Sample),chr.change] <- "n.d."
  }
  
  to.plot.[to.plot.=="early late"] <- "Early"
  to.plot.[to.plot.=="early"] <- "Early"
  to.plot.[to.plot.==" late"] <- "Late"
  to.plot.[to.plot.=="early tetraploid"] <- "n.d."
  to.plot.[to.plot.==" late tetraploid"] <- "n.d."
  to.plot.$`1p`[to.plot.$`1p` %in% c("Early", "Late")] <- "n.d."
  to.plot.$`11q`[to.plot.$`11q` %in% c("Early", "Late")] <- "n.d."
  
  ## equal sizes for boxes
  to.plot.$xmax <- (1:nrow(to.plot.))/nrow(to.plot.)*max.mutation.time.primary/3.3/10^3
  to.plot.$xmin <- ((1:nrow(to.plot.))-1)/nrow(to.plot.)*max.mutation.time.primary/3.3/10^3
  
  plate.subset.list[[length(plate.subset.list) + 1]] <- ggplot(data=to.plot.,
                                                               aes(xmin=xmin, xmax=xmax,
                                                                   ymin=ymin, ymax=ymax,
                                                                   fill=as.character(Telomere.maintenance.mechanism))) +
    geom_rect() +
    geom_rect(data=to.plot., aes(xmin= xmin, xmax=xmax, ymin=ymin+0.27, ymax=ymin+0.48, fill=ManualScore), inherit.aes = F) +
    geom_rect(data=to.plot., aes(xmin= xmin, xmax=xmax, ymin=ymin+0.5, ymax=ymin+0.73, fill=as.character(Rounded.ploidy)), inherit.aes = F) +
    scale_fill_manual(values=c(telomere.colors, clinical.risk.colors, manual.colors, "Subclonal"="purple", "n.d."="grey",
                               ploidy.cols)) +
    scale_x_continuous(limits=c(0, max.mutation.time.primary/3.3/10^3)) +
    theme(legend.position = "bottom") +
    scale_y_continuous(breaks=seq(0.125, 3.125, 0.25),
                       labels=c("TMM", "Stage", "Ploidy", rev(c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2"))))
  
  for(chr.change in rev(c("`17q`", "`1p`", "`1q`", "`7q`", "`2p`", "`11q`", "`17`", "`1`", "`7`", "`2`"))){
    
    index <- which(rev(c("`17q`", "`1p`", "`1q`", "`7q`", "`2p`", "`11q`", "`17`", "`1`", "`7`", "`2`"))==chr.change)
    to.plot.$ymin <- 0.75 + 0.25*(index - 1)+0.2
    to.plot.$ymax <- 0.75 + 0.25*index-0.2
    
    plate.subset.list[[length(plate.subset.list)]] <- plate.subset.list[[length(plate.subset.list)]] +
      geom_rect(data=to.plot., aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill=chr.change), inherit.aes=T)
  }
  
  if(ECA.exists){
    panel="e_ECA_sample_info"
  }else{
    panel="e_no_ECA_sample_info"
  }
  
  addWorksheet(wb, panel)
  to.plot.$Stage <- replace(to.plot.$ManualScore, to.plot.$ManualScore=="LR", "1, 2, 4S")
  to.plot.$Stage <- replace(to.plot.$Stage, to.plot.$Stage=="IR", "3")
  to.plot.$Stage <- replace(to.plot.$Stage, to.plot.$Stage=="HR", "4")
  writeData(wb, panel, to.plot.[,c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2",
                                   "Telomere.maintenance.mechanism", "Stage", "Rounded.ploidy")])
  
}

## are tumors with MNA enriched in the group without ECA

n.total.per.group <- table(c(sample.information.validation[sample.information.validation$Sample.type %in% c("Primary", "Metastasis") &
                                        mutation.time.mrca[rownames(sample.information.validation),]$Mean/3.3/10^3>cutpoint,]$ECA.exists,
                             sample.information.discovery[sample.information.discovery$Sample.type %in% c("Primary", "Metastasis") &
                                                             mutation.time.mrca[rownames(sample.information.discovery),]$Mean/3.3/10^3>cutpoint,]$ECA.exists))


n.mna.per.group <- table(c(sample.information.validation[sample.information.validation$Sample.type %in% c("Primary", "Metastasis") &
                                                         sample.information.validation$Telomere.maintenance.mechanism=="MNA"&
                                                           mutation.time.mrca[rownames(sample.information.validation),]$Mean/3.3/10^3>cutpoint,]$ECA.exists,
                         sample.information.discovery[sample.information.discovery$Sample.type %in% c("Primary", "Metastasis") &
                                                        sample.information.discovery$Telomere.maintenance.mechanism=="MNA"&
                                                     mutation.time.mrca[rownames(sample.information.discovery),]$Mean/3.3/10^3>cutpoint,]$ECA.exists))

test.matrix <- as.matrix(cbind(n.mna.per.group, n.total.per.group - n.mna.per.group))
colnames(test.matrix) <- c("MNA", "no MNA")
fisher.test(test.matrix)

## Plot the cases without ECA separately per ploidy; merge di- and tetraploids
## split them in supplement but not in main figure
pearly.unsplit <- list()
pearly.unsplit.subset.list <- list()

for(ECA.exists in c( F)){
  
  subset=sample.information.validation[ sample.information.validation$Sample.type %in% c("Primary", "Metastasis") &
                                   mutation.time.mrca[rownames(sample.information.validation),]$Mean/3.3/10^3 < cutpoint,,drop=F]
  
  if(nrow(subset)==0){next}
  
  
  subset$Telomere.maintenance.mechanism <- factor(subset$Telomere.maintenance.mechanism,
                                                  levels=c("MNA", "TERT", "ALT", "Multiple", "None"))
  
  subset$MRCAtime <- mutation.time.mrca[rownames(subset),]$Mean
  subset <- subset[order(subset$MRCAtime),]
  
  
  to.plot <- cbind(subset, data.frame(MRCA=mutation.time.mrca[rownames(subset),]$Mean/3.3/10^3,
                                      ECA=mutation.time.eca[rownames(subset),]$Mean/3.3/10^3,
                                      MRCA.upper=mutation.time.mrca[rownames(subset),]$Max/3.3/10^3,
                                      MRCA.lower=mutation.time.mrca[rownames(subset),]$Min/3.3/10^3,
                                      ECA.upper=mutation.time.eca[rownames(subset),]$Max/3.3/10^3,
                                      ECA.lower=mutation.time.eca[rownames(subset),]$Min/3.3/10^3))
  
  
  panel="d"
  
  addWorksheet(wb, panel)
  writeData(wb, panel, to.plot[,c("MRCA", "MRCA.lower", "MRCA.upper")])
  
  
  
  pearly.unsplit[[length(pearly.unsplit)+1]] <- ggplot(data = to.plot[order(to.plot$MRCA),],
                                                       aes(x=MRCA, y=seq(1/length(MRCA),1,length.out = length(MRCA)),
                                                           ymin =  sapply(sort(MRCA), function(x){
                                                             sum(MRCA.upper <= x)
                                                           })/length(MRCA),
                                                           ymax= sapply(sort(MRCA), function(x){
                                                             sum(MRCA.lower <= x)
                                                           })/length(MRCA)
                                                       )) +
    stat_ecdf(col=unname(manual.colors["Late"])) +
    geom_stepribbon(fill=unname(manual.colors["Late"]), alpha=0.5, col=NA)+
    scale_x_continuous(name = "Mutations/Mb",
                       limits=c(0, max.mutation.time.primary/3.3/10^3))+
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(name = "Fraction of tumors") +
    ggtitle(paste("# Cases = ", nrow(subset)))
  
  ## TMM, Clinical subtype & CNVs
  to.plot. <- subset[,c("MRCAtime", "Telomere.maintenance.mechanism", "ManualScore", "Rounded.ploidy")]
  to.plot.$xmin <-c(0, (to.plot.$MRCAtime/3.3/10^3)[-length(to.plot.$MRCAtime)])
  to.plot.$xmax <- to.plot.$MRCAtime/3.3/10^3
  to.plot.$ymin <- 0
  to.plot.$ymax <- 0.23
  
  to.plot.[,c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2")] <- t(matrix.early.late[c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2"),
                                                                                                 rownames(to.plot.)])
  
  ## distinguish cases that are not dateable from cases that are statistically insigniicant from ECA
  for(chr.change in c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2")){
    tmp <- not.dateable[not.dateable$Gain==chr.change,]
    to.plot.[intersect(rownames(to.plot.), tmp$Sample),chr.change] <- "n.d."
  }
  
  to.plot.[to.plot.=="early late"] <- "Early"
  to.plot.[to.plot.=="early"] <- "Early"
  to.plot.[to.plot.==" late"] <- "Late"
  to.plot.[to.plot.=="early tetraploid"] <- "n.d."
  to.plot.[to.plot.==" late tetraploid"] <- "n.d."
  to.plot.$`1p`[to.plot.$`1p` %in% c("Early", "Late")] <- "n.d."
  to.plot.$`11q`[to.plot.$`11q` %in% c("Early", "Late")] <- "n.d."
  
  ## equal sizes for boxes
  to.plot.$xmax <- (1:nrow(to.plot.))/nrow(to.plot.)*max.mutation.time.primary/3.3/10^3
  to.plot.$xmin <- ((1:nrow(to.plot.))-1)/nrow(to.plot.)*max.mutation.time.primary/3.3/10^3
  
  pearly.unsplit.subset.list[[length(pearly.unsplit.subset.list) + 1]] <- ggplot(data=to.plot.,
                                                                                 aes(xmin=xmin, xmax=xmax,
                                                                                     ymin=ymin, ymax=ymax,
                                                                                     fill=as.character(Telomere.maintenance.mechanism))) +
    geom_rect() +
    geom_rect(data=to.plot., aes(xmin= xmin, xmax=xmax, ymin=ymin+0.27, ymax=ymin+0.48, fill=ManualScore), inherit.aes = F) +
    geom_rect(data=to.plot., aes(xmin= xmin, xmax=xmax, ymin=ymin+0.5, ymax=ymin+0.73, fill=as.character(Rounded.ploidy)), inherit.aes = F) +
    scale_fill_manual(values=c(telomere.colors, clinical.risk.colors, manual.colors, "Subclonal"="purple", "n.d."="grey",
                               ploidy.cols)) +
    scale_x_continuous(limits=c(0, max.mutation.time.primary/3.3/10^3)) +
    theme(legend.position = "bottom") +
    scale_y_continuous(breaks=seq(0.125, 3.125, 0.25),
                       labels=c("TMM", "Stage", "Ploidy", rev(c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2"))))
  
  for(chr.change in rev(c("`17q`", "`1p`", "`1q`", "`7q`", "`2p`", "`11q`", "`17`", "`1`", "`7`", "`2`"))){
    
    index <- which(rev(c("`17q`", "`1p`", "`1q`", "`7q`", "`2p`", "`11q`", "`17`", "`1`", "`7`", "`2`"))==chr.change)
    to.plot.$ymin <- 0.75 + 0.25*(index - 1)+0.2
    to.plot.$ymax <- 0.75 + 0.25*index-0.2
    
    pearly.unsplit.subset.list[[length(pearly.unsplit.subset.list)]] <- pearly.unsplit.subset.list[[length(pearly.unsplit.subset.list)]] +
      geom_rect(data=to.plot., aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill=chr.change), inherit.aes=T)
  }
  
  panel="d_sample_info"
  
  addWorksheet(wb, panel)
  to.plot.$Stage <- replace(to.plot.$ManualScore, to.plot.$ManualScore=="LR", "1, 2, 4S")
  to.plot.$Stage <- replace(to.plot.$Stage, to.plot.$Stage=="IR", "3")
  to.plot.$Stage <- replace(to.plot.$Stage, to.plot.$Stage=="HR", "4")
  writeData(wb, panel, to.plot.[,c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2",
                                   "Telomere.maintenance.mechanism", "Stage", "Rounded.ploidy")]) 
  
}


pdf(paste0(panel.directory, "Figure_4_d_e.pdf"), width=9, height=9, useDingbats = F)

figure <- ggarrange(plotlist=c(pearly.unsplit, plate), nrow=3, ncol=3)
annotate_figure(figure, top="Primary tumor / Metastasis")

figure <- ggarrange(plotlist=c(pearly.unsplit.subset.list, plate.subset.list), nrow=3, ncol=3)
annotate_figure(figure, top="Primary tumor / Metastasis")

dev.off()


##########################################################################################################################################
## Figure 4 f,g: Survival curves, overall survival for validation cohort only

addWorksheet(wb, "f")
writeData(wb, "f",  data.frame(Time=EFS.fit.validation$time, n.Risk = EFS.fit.validation$n.risk, 
                                            Category=c(rep(names(EFS.fit.validation$strata)[1], EFS.fit.validation$strata[1]),
                                                       rep(names(EFS.fit.validation$strata)[2], EFS.fit.validation$strata[2])),
                                            Censored=EFS.fit.validation$n.censor, Event = EFS.fit.validation$n.event, Survival=EFS.fit.validation$surv,
                                            LowerCI=EFS.fit.validation$lower, UpperCI=EFS.fit.validation$upper))

chars <- capture.output(EFS.fit.validation.stats)

writeData(wb, sheet = "f", chars, startCol = 10)

addWorksheet(wb, "g")
writeData(wb, "g", data.frame(Time=survival.fit.validation$time, n.Risk = survival.fit.validation$n.risk, 
                              Category=c(rep(names(survival.fit.validation$strata)[1], survival.fit.validation$strata[1]),
                                         rep(names(survival.fit.validation$strata)[2], survival.fit.validation$strata[2])),
                              Censored=survival.fit.validation$n.censor, Event = survival.fit.validation$n.event, Survival=survival.fit.validation$surv,
                              LowerCI=survival.fit.validation$lower, UpperCI=survival.fit.validation$upper))

chars <- capture.output(survival.fit.validation.stats)

writeData(wb, sheet = "g", chars, startCol = 10)


pdf(paste0(panel.directory,"Figure_4f_g.pdf"), useDingbats = F)

ggsurvplot(EFS.fit.validation, data = categorized.by.MRCA.validation, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0,10)) 

ggsurvplot(survival.fit.validation, data = categorized.by.MRCA.validation, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0,10)) 

dev.off()


##########################################################################################################################################
## Figure 4h, i: Survival curves, overall survival, event-free survival

addWorksheet(wb, "h")
writeData(wb, "h", data.frame(Time=joined.EFS.fit$time, n.Risk = joined.EFS.fit$n.risk, 
                               Category=c(rep(names(joined.EFS.fit$strata)[1], joined.EFS.fit$strata[1]),
                                          rep(names(joined.EFS.fit$strata)[2], joined.EFS.fit$strata[2])),
                               Censored=joined.EFS.fit$n.censor, Event = joined.EFS.fit$n.event, Survival=joined.EFS.fit$surv,
                               LowerCI=joined.EFS.fit$lower, UpperCI=joined.EFS.fit$upper))

chars <- capture.output(joined.EFS.fit.stats)

writeData(wb, sheet = "h", chars, startCol = 10)

pdf(paste0(panel.directory,"Figure_4h.pdf"), useDingbats = F)

ggsurvplot(joined.EFS.fit, data = joined.categorized.by.MRCA, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0, 10), break.x.by=5, cex=8, linewidth=0.8) 

dev.off()


addWorksheet(wb, "i")
writeData(wb, "i", data.frame(Time=joined.survival.fit$time, n.Risk = joined.survival.fit$n.risk, 
                               Category=c(rep(names(joined.survival.fit$strata)[1], joined.survival.fit$strata[1]),
                                          rep(names(joined.survival.fit$strata)[2], joined.survival.fit$strata[2])),
                               Censored=joined.survival.fit$n.censor, Event = joined.survival.fit$n.event, Survival=joined.survival.fit$surv,
                               LowerCI=joined.survival.fit$lower, UpperCI=joined.survival.fit$upper))

chars <- capture.output(joined.survival.fit.stats)

writeData(wb, sheet = "i", chars, startCol = 10)

pdf(paste0(panel.directory,"Figure_4i.pdf"), useDingbats = F)
ggsurvplot(joined.survival.fit, data = joined.categorized.by.MRCA, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0, 10), break.x.by=5, cex=8, linewidth=0.8) 
dev.off()



##############################################################################################################################################
## Figure 4j, S4a: Cox regression

pdf(paste0(panel.directory,"Figure_4j.pdf"), useDingbats = F, width = 4, height=4)

ggforest(fit.coxph_MRCA_TMM_Stage_Age_RAS.EFS, data = joined.categorized.by.MRCA, main = "EFS")

dev.off()

chars <- capture.output(summary(fit.coxph_MRCA_TMM_Stage_Age_RAS.EFS))

addWorksheet(wb, "j")
writeData(wb, sheet = "j", chars)


pdf(paste0(panel.directory,"Figure_S4a.pdf"), useDingbats = F, width = 4, height=4)

ggforest(fit.coxph_MRCA_TMM_Stage_Age_RAS.OS, data = joined.categorized.by.MRCA, main="OS")

dev.off()


chars <- capture.output(summary(fit.coxph_MRCA_TMM_Stage_Age_RAS.OS))

addWorksheet(wb.s, "a")
writeData(wb.s, sheet = "a", chars)


##############################################################################################################################################
## Figure S4b/c: Cox regression with RNA classifier


pdf(paste0(panel.directory,"Figure_S4b.pdf"), useDingbats = F, width = 4, height=4)

ggforest(fit.coxph.efs_RNA, data = joined.categorized.by.MRCA[!is.na(joined.categorized.by.MRCA$RNA_classifier),], main = "EFS")

dev.off()

chars <- capture.output(summary(fit.coxph.efs_RNA))

addWorksheet(wb.s, "b")
writeData(wb.s, sheet = "b", chars)


pdf(paste0(panel.directory,"Figure_S4c.pdf"), useDingbats = F, width = 4, height=4)

ggforest(fit.coxph_RNA, data = joined.categorized.by.MRCA[!is.na(joined.categorized.by.MRCA$RNA_classifier),], main="OS")

dev.off()

chars <- capture.output(summary(fit.coxph_RNA))

addWorksheet(wb.s, "c")
writeData(wb.s, sheet = "c", chars)

##########################################################################################################################################

saveWorkbook(wb, file = paste0(panel.directory,"Source_data_Fig.4.xlsx"), overwrite=T)
saveWorkbook(wb.s, file = paste0(panel.directory,"Source_data_Fig.S4.xlsx"), overwrite=T)
