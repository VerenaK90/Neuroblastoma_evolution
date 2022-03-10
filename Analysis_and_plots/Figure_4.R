## Reproduce Fig. 4
##############################################################################################################################################
## load settings and libraries

load(paste0(rdata.directory, "MRCA_timing.RData"))

## source data:
wb <- createWorkbook()

## store figure panels
panel.directory <- paste0(output.directory, "Figure4/")

if(!dir.exists(panel.directory)){
  dir.create(panel.directory)
}


##########################################################################################################################################
## Figure 4a, b: Survival curves, overall survival for discovery cohort only

## EFS
addWorksheet(wb.s, "a")
writeData(wb.s, "a",  data.frame(Time=EFS.fit$time, n.Risk = EFS.fit$n.risk, 
                                 Category=c(rep(names(EFS.fit$strata)[1], EFS.fit$strata[1]),
                                            rep(names(EFS.fit$strata)[2], EFS.fit$strata[2])),
                                 Censored=EFS.fit$n.censor, Event = EFS.fit$n.event, Survival=EFS.fit$surv,
                                 LowerCI=EFS.fit$lower, UpperCI=EFS.fit$upper))

## OS
addWorksheet(wb, "b")
writeData(wb, "b", data.frame(Time=survival.fit$time, n.Risk = survival.fit$n.risk, 
                              Category=c(rep(names(survival.fit$strata)[1], survival.fit$strata[1]),
                                         rep(names(survival.fit$strata)[2], survival.fit$strata[2])),
                              Censored=survival.fit$n.censor, Event = survival.fit$n.event, Survival=survival.fit$surv,
                              LowerCI=survival.fit$lower, UpperCI=survival.fit$upper))


pdf(paste0(panel.directory,"Figure_4a_b.pdf"), useDingbats = F)

ggsurvplot(EFS.fit, data = categorized.by.MRCA, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0,10)) 

ggsurvplot(survival.fit, data = categorized.by.MRCA, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0,10)) 

dev.off()


##########################################################################################################################################
## Figure 4c, Mutation density at MRCA in the validation cohort


pdf(paste0(panel.directory, "Figure_4c.pdf"), width=3.5, height=3.5, useDingbats = F)

to.plot <- data.frame(MRCA=mutation.time.mrca[primary.tumors.30x])
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

source(paste0(custom.script.directory, "Oncoprint_additional_samples.R"))

matrix.early.late <- matrix.early.late.30x
## if ECA was not uniquely identified, take the earliest time point
mutation.time.eca[names(earliest.mutation.time)] <- earliest.mutation.time
mutation.time.eca.lower[names(earliest.mutation.time)] <- earliest.mutation.time.lower
mutation.time.eca.upper[names(earliest.mutation.time)] <- earliest.mutation.time.upper
####

sample.information.30x$Ploidy <- sample.information.30x$Rounded.ploidy

max.mutation.time.primary <- max(mutation.time.mrca[rownames(sample.information.30x[sample.information.30x$Location %in% c("Primary", "Relapse"),])],
                                 na.rm=T)

sample.information.30x$Telomere.maintenance.mechanism <- factor(sample.information.30x$Telomere.maintenance.mechanism,
                                                                levels=c("MNA", "TERT", "ALT", "Multiple", "None"))

sample.information.30x$ECA.exists <- as.logical(sample.information.30x$ECA.exists)

##### Primary tumors, with and without treatment, metastases
## Plot early and late MRCA separately; further, stratify late cases w.r.t. number of events

cutpoint <- 0.05

plate <- list()
plate.subset.list <- list()

for(ECA.exists in c(T, F)){
  
  subset=sample.information.30x[ sample.information.30x$Location %in% c("Primary", "Metastasis") &
                                   sample.information.30x$ECA.exists==ECA.exists &
                                   mutation.time.mrca[rownames(sample.information.30x)]/3.3/10^3 >=cutpoint,,drop=F]
  
  if(nrow(subset)==0){next}
  
  
  subset$Telomere.maintenance.mechanism <- factor(subset$Telomere.maintenance.mechanism,
                                                  levels=c("MNA", "TERT", "ALT", "Multiple", "None"))
  
  subset$MRCAtime <- mutation.time.mrca[rownames(subset)]
  subset <- subset[order(subset$MRCAtime),]
  
  
  to.plot <- cbind(subset, data.frame(MRCA=mutation.time.mrca[rownames(subset)]/3.3/10^3,
                                      ECA=mutation.time.eca[rownames(subset)]/3.3/10^3,
                                      MRCA.upper=mutation.time.mrca.upper[rownames(subset)]/3.3/10^3,
                                      MRCA.lower=mutation.time.mrca.lower[rownames(subset)]/3.3/10^3,
                                      ECA.upper=mutation.time.eca.upper[rownames(subset)]/3.3/10^3,
                                      ECA.lower=mutation.time.eca.lower[rownames(subset)]/3.3/10^3))
  
  
  
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
  to.plot. <- subset[,c("MRCAtime", "Telomere.maintenance.mechanism", "ManualScore", "Ploidy")]
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
  
  plate.subset.list[[length(plate.subset.list) + 1]] <- ggplot(data=to.plot.,
                                                               aes(xmin=xmin, xmax=xmax,
                                                                   ymin=ymin, ymax=ymax,
                                                                   fill=as.character(Telomere.maintenance.mechanism))) +
    geom_rect() +
    geom_rect(data=to.plot., aes(xmin= xmin, xmax=xmax, ymin=ymin+0.27, ymax=ymin+0.48, fill=ManualScore), inherit.aes = F) +
    geom_rect(data=to.plot., aes(xmin= xmin, xmax=xmax, ymin=ymin+0.5, ymax=ymin+0.73, fill=as.character(Ploidy)), inherit.aes = F) +
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
                                   "Telomere.maintenance.mechanism", "Stage", "Ploidy")])
  
}


## Plot the cases without ECA separately per ploidy; merge di- and tetraploids
## split them in supplement but not in main figure
pearly.unsplit <- list()
pearly.unsplit.subset.list <- list()

for(ECA.exists in c( F)){
  
  subset=sample.information.30x[ sample.information.30x$Location %in% c("Primary", "Metastasis") &
                                   mutation.time.mrca[rownames(sample.information.30x)]/3.3/10^3 < cutpoint,,drop=F]
  
  if(nrow(subset)==0){next}
  
  
  subset$Telomere.maintenance.mechanism <- factor(subset$Telomere.maintenance.mechanism,
                                                  levels=c("MNA", "TERT", "ALT", "Multiple", "None"))
  
  subset$MRCAtime <- mutation.time.mrca[rownames(subset)]
  subset <- subset[order(subset$MRCAtime),]
  
  
  to.plot <- cbind(subset, data.frame(MRCA=mutation.time.mrca[rownames(subset)]/3.3/10^3,
                                      ECA=mutation.time.eca[rownames(subset)]/3.3/10^3,
                                      MRCA.upper=mutation.time.mrca.upper[rownames(subset)]/3.3/10^3,
                                      MRCA.lower=mutation.time.mrca.lower[rownames(subset)]/3.3/10^3,
                                      ECA.upper=mutation.time.eca.upper[rownames(subset)]/3.3/10^3,
                                      ECA.lower=mutation.time.eca.lower[rownames(subset)]/3.3/10^3))
  
  
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
  to.plot. <- subset[,c("MRCAtime", "Telomere.maintenance.mechanism", "ManualScore", "Ploidy")]
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
    geom_rect(data=to.plot., aes(xmin= xmin, xmax=xmax, ymin=ymin+0.5, ymax=ymin+0.73, fill=as.character(Ploidy)), inherit.aes = F) +
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
                                   "Telomere.maintenance.mechanism", "Stage", "Ploidy")]) 
  
}


pdf(paste0(panel.directory, "Figure_4_d_e.pdf"), width=9, height=9, useDingbats = F)

figure <- ggarrange(plotlist=c(pearly.unsplit, plate), nrow=3, ncol=3)
annotate_figure(figure, top="Primary tumor / Metastasis")

figure <- ggarrange(plotlist=c(pearly.unsplit.subset.list, plate.subset.list), nrow=3, ncol=3)
annotate_figure(figure, top="Primary tumor / Metastasis")

dev.off()

##########################################################################################################################################
## Figure 4h, i: Survival curves, overall survival, event-free survival

addWorksheet(wb, "4h")
writeData(wb, "4h", data.frame(Time=joined.EFS.fit$time, n.Risk = joined.EFS.fit$n.risk, 
                               Category=c(rep(names(joined.EFS.fit$strata)[1], joined.EFS.fit$strata[1]),
                                          rep(names(joined.EFS.fit$strata)[2], joined.EFS.fit$strata[2])),
                               Censored=joined.EFS.fit$n.censor, Event = joined.EFS.fit$n.event, Survival=joined.EFS.fit$surv,
                               LowerCI=joined.EFS.fit$lower, UpperCI=joined.EFS.fit$upper))

pdf(paste0(panel.directory,"Figure_4h.pdf"), useDingbats = F)

ggsurvplot(joined.EFS.fit, data = joined.categorized.by.MRCA, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0, 10)) 

dev.off()


addWorksheet(wb, "4i")
writeData(wb, "4i", data.frame(Time=joined.survival.fit$time, n.Risk = joined.survival.fit$n.risk, 
                               Category=c(rep(names(joined.survival.fit$strata)[1], joined.survival.fit$strata[1]),
                                          rep(names(joined.survival.fit$strata)[2], joined.survival.fit$strata[2])),
                               Censored=joined.survival.fit$n.censor, Event = joined.survival.fit$n.event, Survival=joined.survival.fit$surv,
                               LowerCI=joined.survival.fit$lower, UpperCI=joined.survival.fit$upper))

pdf(paste0(panel.directory,"Figure_4i.pdf"), useDingbats = F)
ggsurvplot(joined.survival.fit, data = joined.categorized.by.MRCA, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0, 10)) 
dev.off()




##############################################################################################################################################
## Figure 4j: Cox regression


pdf(paste0(panel.directory,"Figure_4j.pdf"), useDingbats = F, width = 4, height=4)

ggforest(fit.coxph.efs, data = joined.categorized.by.MRCA, main = "EFS")

dev.off()
##########################################################################################################################################

saveWorkbook(wb, file = paste0(panel.directory,"Source_data_Fig.4.xlsx"), overwrite=T)
