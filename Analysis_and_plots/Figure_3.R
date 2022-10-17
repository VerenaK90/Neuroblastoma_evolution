## Reproduce Fig. 3
##############################################################################################################################################
## load settings and libraries
source(".Settings.R")

#source(paste0(custom.script.directory, "Oncoprint.R"))
#source(paste0(custom.script.directory, "Mutation_densities_at_chromosomal_gains.R"))
load(paste0(rdata.directory, "MRCA_timing.RData"))
#load(paste0(rdata.directory, "Estimated_mutation_rate_per_day.RData"))
load(paste0(rdata.directory, "Sig_colors.RData"))
load(paste0(rdata.directory, "Clock_like_SBS.RData"))
load(paste0(rdata.directory, "Clonal_mutations_different_CNs.RData"))

## source data:
wb <- createWorkbook()
## add some panels to Supplementary Fig. 2, some to 3
wb.s2 <- loadWorkbook(paste0(paste0(output.directory, "Figure2/"), "Source_data_Fig.S2.xlsx"))
wb.s3 <- loadWorkbook(paste0(paste0(output.directory, "Figure2/"), "Source_data_Fig.S3.xlsx"))

## store figure panels
panel.directory <- paste0(output.directory, "Figure3/")

if(!dir.exists(panel.directory)){
  dir.create(panel.directory)
}

## MRCA cutpoint
cutpoint <- 0.05

##############################################################################################################################################
## Figure 3a: Mutation density at MRCAs in the discovery set
## color by cutpoint


pdf(paste0(panel.directory, "Figure_3a.pdf"), width=3.5, height=3.5, useDingbats = F)

to.plot <- data.frame(MRCA=mutation.time.mrca[primary.tumors.discovery,]$Mean)
to.plot$MRCA.time <- ifelse(to.plot$MRCA/3.3/10^3 < 0.05, "Early", "Late")

addWorksheet(wb, "a")
writeData(wb, "a", data.frame(MRCA=to.plot$MRCA/3.3/10^3,
                              MRCA.time=to.plot$MRCA.time))


ggplot(to.plot, aes(x=MRCA/3.3/10^3, fill=MRCA.time)) + geom_histogram(binwidth=0.02) +
  scale_fill_manual(values=c(Early="dodgerblue", Late= "dodgerblue4")) +
  geom_vline(xintercept = cutpoint, linetype=2) + 
  scale_x_continuous(name="#SSNVs/Mb at MRCA") + scale_y_continuous(name="# Tumors") +
  ggtitle("Discovery set")

dev.off()

##############################################################################################################################################
## Figure S2d/e: Is the difference in mutation density due to a different mutational process? 

## clonal SNVs only
clonal.sigs.discovery <- read.delim(paste0(signature.directory, "/discovery/MMSig_output_discovery_C.tsv"))
clonal.sigs.discovery$Clock <- clonal.sigs.discovery$SBS1 + clonal.sigs.discovery$SBS5 + clonal.sigs.discovery$SBS40
clonal.sigs.discovery$Sample <- rownames(clonal.sigs.discovery)
clonal.sigs.discovery$MRCA.time <- replace(rep("Late", nrow(clonal.sigs.discovery)), mutation.time.mrca[clonal.sigs.discovery$Sample ,]$Mean/3.3/10^3 < cutpoint, "Early")
clonal.sigs.discovery$Diagnosis <- replace(rep("Primary", nrow(clonal.sigs.discovery)), !sample.information.discovery[clonal.sigs.discovery$Sample,]$Sample.type %in% c("Primary", "Metastasis"), "Relapse")

to.plot <- melt(clonal.sigs.discovery, value.name = "Exposure", id.vars = c("Sample", "MRCA.time", "Diagnosis"), 
                variable.name = "Signature")
to.plot$Signature <- factor(to.plot$Signature, levels = c("Clock", setdiff(unique(to.plot$Signature), "Clock")))

p1 <- ggplot(to.plot[!to.plot$Signature %in% c("SBS1", "SBS5", "SBS40", "mutations"),], aes(x=Signature, y=Exposure)) + 
  geom_bar(stat="summary") + geom_errorbar(stat = "summary") + facet_wrap(~MRCA.time + Diagnosis)

to.store <- to.plot[!to.plot$Signature %in% c("SBS1", "SBS5", "SBS40", "mutations"),]
to.store$Sample <- sapply(to.store$Sample, function(x){
  if(x %in% tumors.discovery){
    sample.information.discovery[x,]$Evolution_paper_Id
  }else{
    sample.information.validation[x,]$Evolution_paper_ID
  }
})

addWorksheet(wb.s2, "d")
writeData(wb.s2, "d", to.store)

## subclonal SNVs only
subclonal.sigs.discovery <- read.delim(paste0(signature.directory, "/discovery/MMSig_output_discovery_SC.tsv"))
subclonal.sigs.discovery$Clock <- subclonal.sigs.discovery$SBS1 + subclonal.sigs.discovery$SBS5 + subclonal.sigs.discovery$SBS40
subclonal.sigs.discovery$Sample <- rownames(subclonal.sigs.discovery)
subclonal.sigs.discovery$MRCA.time <- replace(rep("Late", nrow(subclonal.sigs.discovery)), mutation.time.mrca[subclonal.sigs.discovery$Sample ,]$Mean/3.3/10^3 < cutpoint, "Early")
subclonal.sigs.discovery$Diagnosis <- replace(rep("Primary", nrow(subclonal.sigs.discovery)), !sample.information.discovery[subclonal.sigs.discovery$Sample,]$Sample.type %in% c("Primary", "Metastasis"), "Relapse")

to.plot <- melt(subclonal.sigs.discovery, value.name = "Exposure", id.vars = c("Sample", "MRCA.time", "Diagnosis"), 
                variable.name = "Signature")
to.plot$Signature <- factor(to.plot$Signature, levels = c("Clock", setdiff(unique(to.plot$Signature), "Clock")))

p2 <- ggplot(to.plot[!to.plot$Signature %in% c("SBS1", "SBS5", "SBS40", "mutations"),], aes(x=Signature, y=Exposure)) + 
  geom_bar(stat="summary") + geom_errorbar(stat = "summary") + facet_wrap(~MRCA.time + Diagnosis)


pdf(paste0(panel.directory, "Figure_S2d_e.pdf"), width=3.5, height=3.5, useDingbats = F)

p1 + ggtitle("Clonal")

p2 + ggtitle("Subclonal")

dev.off()

to.store <- to.plot[!to.plot$Signature %in% c("SBS1", "SBS5", "SBS40", "mutations"),]
to.store$Sample <- sapply(to.store$Sample, function(x){
  if(x %in% tumors.discovery){
    sample.information.discovery[x,]$Evolution_paper_Id
  }else{
    sample.information.validation[x,]$Evolution_paper_ID
  }
})

addWorksheet(wb.s2, "e")
writeData(wb.s2, "e", to.store)

##############################################################################################################################################
## Figure 3b,c: density distributions for MRCA and ECA - discovery set
source(paste0(custom.script.directory, "Driver_mutations.R"))
source(paste0(custom.script.directory, "Oncoprint.R"))

## if ECA was not uniquely identified, take the earliest time point
mutation.time.eca[earliest.mutation.time$Sample,]$Mean <- earliest.mutation.time$Mean
mutation.time.eca[earliest.mutation.time$Sample,]$Min <- earliest.mutation.time$Min
mutation.time.eca[earliest.mutation.time$Sample,]$Max <- earliest.mutation.time$Max
####

colnames(sample.information.discovery)[which(colnames(sample.information.discovery)=="ECA")] <- "ECA.exists"


max.mutation.time.primary <- max(mutation.time.mrca[rownames(sample.information.discovery[sample.information.discovery$Sample.type %in% c("Primary", "Relapse"),]),]$Max)

sample.information.discovery$Telomere.maintenance.mechanism <- factor(sample.information.discovery$Telomere.maintenance.mechanism,
                                                                levels=c("MNA", "TERT", "ALT", "Multiple", "None"))

sample.information.discovery$ECA.exists <- as.logical(sample.information.discovery$ECA.exists)

##### Primary tumors, with and without treatment, metastases
## Plot early and late MRCA separately; further, stratify late cases w.r.t. number of events


plate <- list()
plate.subset.list <- list()

for(ECA.exists in c(T, F)){
  
  subset=sample.information.discovery[ sample.information.discovery$Sample.type %in% c("Primary", "Metastasis") &
                                   sample.information.discovery$ECA.exists==ECA.exists &
                                   mutation.time.mrca[rownames(sample.information.discovery),]$Mean/3.3/10^3 >=cutpoint,,drop=F]
  
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
    panel="c_ECA"
    mean(to.plot$MRCA - to.plot$ECA)
    ## 0.24
    sd(to.plot$MRCA - to.plot$ECA)
    ## 0.1
  }else{
    panel="c_no_ECA"
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
    panel="c_ECA_sample_info"
  }else{
    panel="c_no_ECA_sample_info"
  }
  
  addWorksheet(wb, panel)
  to.plot.$Stage <- replace(to.plot.$ManualScore, to.plot.$ManualScore=="LR", "1, 2, 4S")
  to.plot.$Stage <- replace(to.plot.$Stage, to.plot.$Stage=="IR", "3")
  to.plot.$Stage <- replace(to.plot.$Stage, to.plot.$Stage=="HR", "4")
  writeData(wb, panel, to.plot.[,c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2",
                                   "Telomere.maintenance.mechanism", "Stage", "Rounded.ploidy")])
  
}


## Plot the cases without ECA separately per ploidy; merge di- and tetraploids
## split them in supplement but not in main figure
pearly.unsplit <- list()
pearly.unsplit.subset.list <- list()

for(ECA.exists in c( F)){
  
  subset=sample.information.discovery[ sample.information.discovery$Sample.type %in% c("Primary", "Metastasis") &
                                   mutation.time.mrca[rownames(sample.information.discovery),]$Mean/3.3/10^3 < cutpoint,,drop=F]
  
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
  
  
  panel="b"
  
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
  
  panel="b_sample_info"
  
  addWorksheet(wb, panel)
  to.plot.$Stage <- replace(to.plot.$ManualScore, to.plot.$ManualScore=="LR", "1, 2, 4S")
  to.plot.$Stage <- replace(to.plot.$Stage, to.plot.$Stage=="IR", "3")
  to.plot.$Stage <- replace(to.plot.$Stage, to.plot.$Stage=="HR", "4")
  writeData(wb, panel, to.plot.[,c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2",
                                   "Telomere.maintenance.mechanism", "Stage", "Rounded.ploidy")]) 
  
}


pdf(paste0(panel.directory, "Figure_3_b_c.pdf"), width=9, height=9, useDingbats = F)

figure <- ggarrange(plotlist=c(pearly.unsplit, plate), nrow=3, ncol=3)
annotate_figure(figure, top="Primary tumor / Metastasis")

figure <- ggarrange(plotlist=c(pearly.unsplit.subset.list, plate.subset.list), nrow=3, ncol=3)
annotate_figure(figure, top="Primary tumor / Metastasis")

dev.off()


##########################################################################################################################################
## Figure S33: Compare MRCA between cases with and without ECA among late-MRCA tumors

pdf(paste0(panel.directory,"Figure_S3e.pdf"), useDingbats = F)

to.plot <- data.frame(MRCA=mutation.time.mrca[primary.tumors.discovery,]$Mean/3.3/10^3,
                      ECA.exists=ifelse(is.na(mutation.time.eca[primary.tumors.discovery,]$Mean) | 
                                          mutation.time.eca[primary.tumors.discovery,]$Mean==mutation.time.mrca[primary.tumors.discovery,]$Mean, F, T),
                      MRCA.time=ifelse(mutation.time.mrca[primary.tumors.discovery,]$Mean/3.3/10^3<cutpoint, "early", "late"))

addWorksheet(wb.s3, "e")
writeData(wb.s3, "e", to.plot)

ggplot(to.plot[to.plot$MRCA.time=="late" ,], 
       aes(x=ECA.exists, y=MRCA, col=ECA.exists)) + geom_boxplot() + geom_beeswarm() +
  scale_color_manual(values=c("TRUE" = unname(manual.colors["Early"]), "FALSE" = unname(manual.colors["Late"]))) +
  geom_boxplot(data=to.plot[to.plot$MRCA.time=="early",], col="dodgerblue", aes(x="Early MRCA"))+
  geom_beeswarm(data=to.plot[to.plot$MRCA.time=="early",], col="dodgerblue", aes(x="Early MRCA"))

dev.off()


##########################################################################################################################################
#### Figure S2i: # of events for triploidization

mutation.time.eca$Mean[is.na(mutation.time.eca$Mean)] <- mutation.time.mrca$Mean[is.na(mutation.time.eca$Mean)]

## For each ploidy: can we assign most polyploidization events to a single event? 
numbers.of.3n.chromosomes.at.single.event <- c()
numbers.of.3n.chromosomes <- c()
not.conforming.chromosomes <- list()

for(i in triploid.tumors.discovery){
  
  ## take only dominating copy number per chromosome
  
  dominating.copy.number.per.chromosome <- paste("chr", 1:22, 
                                                 rownames(clonal.mutations.all.tumors[[i]]$segment.length.matrix[-c(4,7),])[
                                                   apply(clonal.mutations.all.tumors[[i]]$segment.length.matrix[-c(4,7),], 2, which.max)],
                                                 sep="_")
    
  
  if(mutation.time.eca[i,]$Mean==mutation.time.mrca[i,]$Mean){
    
    ## if these have a higher mutation count than the mrca estimate count them to mrca
    
    
    if(length(mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca)>0){
      tmp <- c(mrca.eca[[i]]$gains.at.mrca, 
               mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca[sapply(mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca, function(x){
                 mutation.time[[i]][ mutation.time[[i]]$Segment==x,]$MRCA >= mutation.time.mrca[i,]$Mean
               })])
      
      tmp <- sapply(intersect(tmp, dominating.copy.number.per.chromosome), function(x){strsplit(x, split="_")[[1]][3]})
      
      tmp.2 <- mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca[sapply(mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca, function(x){
        mutation.time[[i]][ mutation.time[[i]]$Segment==x,]$Mean < mutation.time.mrca[i,]$Mean
      })]
      tmp.2 <- intersect(tmp.2, dominating.copy.number.per.chromosome)
      tmp.2 <- sapply(intersect(tmp.2, dominating.copy.number.per.chromosome), function(x){strsplit(x, split="_")[[1]][3]})
    }else{
      tmp <- mrca.eca[[i]]$gains.at.mrca
      
      tmp <- sapply(intersect(tmp, dominating.copy.number.per.chromosome), function(x){strsplit(x, split="_")[[1]][3]})
      
      tmp.2 <- c()
    }
    
    numbers.of.3n.chromosomes.at.single.event[i] <- sum(tmp %in% c("trisomic", "tetrasomic"))
    numbers.of.3n.chromosomes[i] <- (sum(tmp.2 %in% c("trisomic", "tetrasomic"))+sum(tmp %in% c("trisomic", "tetrasomic")))
    
    tmp.2 <- sapply(names(tmp.2[tmp.2 %in% c("trisomic", "tetrasomic")]), function(x){strsplit(x, split="_")[[1]][2]})
    not.conforming.chromosomes[[i]] <- tmp.2
    
  }else{
    tmp <- c(mrca.eca[[i]]$gains.uniquely.mapped.to.eca, mrca.eca[[i]]$gains.at.mrca.conforming.eca)
    tmp <- sapply(intersect(tmp, dominating.copy.number.per.chromosome), function(x){strsplit(x, split="_")[[1]][3]})
    
    if(length(mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca)>0){
      tmp.2 <- c(setdiff(c(mrca.eca[[i]]$gains.at.mrca, mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca[sapply(mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca, function(x){
        mutation.time[[i]][ as.character(mutation.time[[i]]$Segment)==x,]$Mean >= mutation.time.mrca[i,]$Mean
      })]), mrca.eca[[i]]$gains.at.mrca.conforming.eca), mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca[sapply(mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca, function(x){
        mutation.time[[i]][ as.character(mutation.time[[i]]$Segment)==x,]$Mean < mutation.time.mrca[i,]$Mean
      })])
    }else{
      tmp.2 <- setdiff(mrca.eca[[i]]$gains.at.mrca, mrca.eca[[i]]$gains.at.mrca.conforming.eca)
    }


    tmp.2 <- intersect(tmp.2, dominating.copy.number.per.chromosome)
    tmp.2 <- sapply(intersect(tmp.2, dominating.copy.number.per.chromosome), function(x){strsplit(x, split="_")[[1]][3]})
    
    numbers.of.3n.chromosomes.at.single.event[i] <- sum(tmp %in% c("trisomic", "tetrasomic"))
    numbers.of.3n.chromosomes[i] <- (sum(tmp.2 %in% c("trisomic", "tetrasomic"))+sum(tmp %in% c("trisomic", "tetrasomic")))
    
    tmp.2 <- sapply(names(tmp.2[tmp.2 %in% c("trisomic", "tetrasomic")]), function(x){strsplit(x, split="_")[[1]][2]})
    not.conforming.chromosomes[[i]] <- tmp.2
  }
}


numbers.of.4n.chromosomes.at.single.event <- c()
numbers.of.4n.chromosomes <- c()
not.conforming.chromosomes <- list()

for(i in tetraploid.tumors.discovery){
  
  ## take only dominating copy number per chromosome
  
  dominating.copy.number.per.chromosome <- paste("chr", 1:22, 
                                                 rownames(clonal.mutations.all.tumors[[i]]$segment.length.matrix[-c(4,7),])[
                                                   apply(clonal.mutations.all.tumors[[i]]$segment.length.matrix[-c(4,7),], 2, which.max)],
                                                 sep="_")
  
  if(mutation.time.eca[i,]$Mean==mutation.time.mrca[i,]$Mean){
    
    ## if these have a higher mutation count than the mrca estimate count them to mrca
    if(length(mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca)>0){
      tmp <- c(mrca.eca[[i]]$gains.at.mrca, 
               mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca[sapply(mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca, function(x){
                 mutation.time[[i]][ mutation.time[[i]]$Segment==x,]$Mean >= mutation.time.mrca[i,]$Mean
               })])
      
      tmp <- sapply(intersect(tmp, dominating.copy.number.per.chromosome), function(x){strsplit(x, split="_")[[1]][3]})
      
      
      
      tmp.2 <- mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca[sapply(mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca, function(x){
        mutation.time[[i]][ mutation.time[[i]]$Segment==x,]$Mean < mutation.time.mrca[i,]$Mean
      })]
      
      
      tmp.2 <- intersect(tmp.2, dominating.copy.number.per.chromosome)
      tmp.2 <- sapply(intersect(tmp.2, dominating.copy.number.per.chromosome), function(x){strsplit(x, split="_")[[1]][3]})
      
    }else{
      tmp <- mrca.eca[[i]]$gains.at.mrca
      
      tmp <- sapply(intersect(tmp, dominating.copy.number.per.chromosome), function(x){strsplit(x, split="_")[[1]][3]})
      
      tmp.2 <- c()
      
    }
      
    
    numbers.of.4n.chromosomes.at.single.event[i] <- sum(tmp %in% c("trisomic", "tetrasomic"))
    numbers.of.4n.chromosomes[i] <- (sum(tmp.2 %in% c("trisomic", "tetrasomic"))+sum(tmp %in% c("trisomic", "tetrasomic")))
    
    tmp.2 <- sapply(names(tmp.2[tmp.2 %in% c("trisomic", "tetrasomic")]), function(x){strsplit(x, split="_")[[1]][2]})
    not.conforming.chromosomes[[i]] <- tmp.2
    
  }else{
    tmp <- c(mrca.eca[[i]]$gains.uniquely.mapped.to.eca, mrca.eca[[i]]$gains.at.mrca.conforming.eca)
    tmp <- sapply(intersect(tmp, dominating.copy.number.per.chromosome), function(x){strsplit(x, split="_")[[1]][3]})
    
    if(length(mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca)>0){
      tmp.2 <- c(setdiff(c(mrca.eca[[i]]$gains.at.mrca, mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca[sapply(mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca, function(x){
        mutation.time[[i]][ as.character(mutation.time[[i]]$Segment)==x,]$Mean >= mutation.time.mrca[i,]$Mean
      })]), mrca.eca[[i]]$gains.at.mrca.conforming.eca), mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca[sapply(mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca, function(x){
        mutation.time[[i]][ as.character(mutation.time[[i]]$Segment)==x,]$Mean < mutation.time.mrca[i,]$Mean
      })])
    }else{
      tmp.2 <- setdiff(mrca.eca[[i]]$gains.at.mrca, mrca.eca[[i]]$gains.at.mrca.conforming.eca)
    }
    
    
    tmp.2 <- intersect(tmp.2, dominating.copy.number.per.chromosome)
    tmp.2 <- sapply(intersect(tmp.2, dominating.copy.number.per.chromosome), function(x){strsplit(x, split="_")[[1]][3]})
    
    numbers.of.4n.chromosomes.at.single.event[i] <- sum(tmp %in% c("trisomic", "tetrasomic"))
    numbers.of.4n.chromosomes[i] <- (sum(tmp.2 %in% c("trisomic", "tetrasomic"))+sum(tmp %in% c("trisomic", "tetrasomic")))
    
    tmp.2 <- sapply(names(tmp.2[tmp.2 %in% c("trisomic", "tetrasomic")]), function(x){strsplit(x, split="_")[[1]][2]})
    not.conforming.chromosomes[[i]] <- tmp.2
  }
}


numbers.of.2n.chromosomes.at.mrca <- c()
numbers.of.2n.chromosomes <- c()
not.conforming.chromosomes <- list()

for(i in diploid.tumors.discovery){
  
  ## take only dominating copy number per chromosome
  
  dominating.copy.number.per.chromosome <- paste("chr", 1:22, 
                                                 rownames(clonal.mutations.all.tumors[[i]]$segment.length.matrix[-c(4,7),])[
                                                   apply(clonal.mutations.all.tumors[[i]]$segment.length.matrix[-c(4,7),], 2, which.max)],
                                                 sep="_")
  
  ## among them select the disomic chromosomes
  
  disomic.chromosomes <- dominating.copy.number.per.chromosome[sapply(dominating.copy.number.per.chromosome, function(x){
    x <- strsplit(x, split="_")[[1]]
    if(x[3]=="disomic"){
      return(T)
    }else{
      return(F)
    }
  })]
  
  tmp <- setdiff(disomic.chromosomes, mrca.eca[[i]]$monosomic.states.not.matching.mrca)
  
  if(length(mrca.eca[[i]]$monosomic.states.not.matching.mrca)==0){
    tmp.2 <- mrca.eca[[i]]$monosomic.states.not.matching.mrca
  }else{
    tmp.2 <- mrca.eca[[i]]$monosomic.states.not.matching.mrca[sapply(mrca.eca[[i]]$monosomic.states.not.matching.mrca, function(x){
      x <- strsplit(x, split="_")[[1]][3]
      if(x=="disomic"){
        return(T)
      }else{
        return(F)
      }
    })]
  }
  
  numbers.of.2n.chromosomes.at.mrca[i] <- length(tmp)
  numbers.of.2n.chromosomes[i] <- length(tmp) + length(tmp.2)
  
  not.conforming.chromosomes[[i]] <- tmp.2
  
}

pdf(paste0(panel.directory,"Figure_S3f.pdf"), width=2, height=2.2, useDingbats = F)

to.plot <- rbind(data.frame(Chromosomes=numbers.of.2n.chromosomes.at.mrca/numbers.of.2n.chromosomes,
                            Tumor=diploid.tumors.discovery,
                            Telomere.type=telomere.classification.discovery[diploid.tumors.discovery],
                            Time = sample.information.discovery[diploid.tumors.discovery,]$Sample.type,
                            Ploidy=2),
                 data.frame(Chromosomes=numbers.of.3n.chromosomes.at.single.event/numbers.of.3n.chromosomes,
                            Tumor=triploid.tumors.discovery,
                            Telomere.type=telomere.classification.discovery[triploid.tumors.discovery],
                            Time = sample.information.discovery[triploid.tumors.discovery,]$Sample.type,
                            Ploidy=3),
                 data.frame(Chromosomes=numbers.of.4n.chromosomes.at.single.event/numbers.of.4n.chromosomes,
                            Tumor=tetraploid.tumors.discovery,
                            Telomere.type=telomere.classification.discovery[tetraploid.tumors.discovery],
                            Time = sample.information.discovery[tetraploid.tumors.discovery,]$Sample.type,
                            Ploidy=4))


panel="f"
addWorksheet(wb.s3, panel)
writeData(wb.s3, panel, to.plot)


ggplot(to.plot, aes(x=Ploidy, y=Chromosomes,group=Ploidy)) + geom_boxplot(width=0.5) + geom_beeswarm() + scale_y_continuous(limits=c(0,1.01))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


dev.off()

##########################################################################################################################################
## Fig. 3d: compare mutation density at ECA of late-MRCA group with MRCA of early-MRCAgroup

pdf(paste0(panel.directory,"Figure_3d.pdf"), useDingbats = F, width=2, height=2)

to.plot <- data.frame(MRCA=mutation.time.mrca[primary.tumors.discovery,]$Mean/3.3/10^3,
                      ECA=mutation.time.eca[primary.tumors.discovery,]$Mean/3.3/10^3,
                      ECA.exists=ifelse(is.na(mutation.time.eca[primary.tumors.discovery,]$Mean) | mutation.time.eca[primary.tumors.discovery,]$Mean!=mutation.time.mrca[primary.tumors.discovery,]$Mean, T, F),
                      MRCA.time=ifelse(mutation.time.mrca[primary.tumors.discovery,]$Mean/3.3/10^3<cutpoint, "early", "late"))

to.plot <- to.plot[to.plot$ECA.exists==T | to.plot$MRCA.time=="early",]

addWorksheet(wb, "d")
writeData(wb, "d", to.plot)

ggplot(to.plot[to.plot$MRCA.time=="late" ,], 
       aes(x=ECA.exists, y=ECA)) + geom_boxplot(col="dodgerblue4") + geom_beeswarm(col="dodgerblue4") +
  geom_boxplot(data=to.plot[to.plot$MRCA.time=="early",], col="dodgerblue", aes(x="Early MRCA", y=MRCA))+
  geom_beeswarm(data=to.plot[to.plot$MRCA.time=="early",], col="dodgerblue", aes(x="Early MRCA", y=MRCA)) +
  scale_y_continuous(name="#SSNVs/Mb")

chars <- capture.output(wilcox.test(to.plot[to.plot$MRCA.time=="late",]$ECA, 
            to.plot[to.plot$MRCA.time=="early",]$MRCA, conf.int = T))


writeData(wb, sheet = "d", chars, startCol = ncol(to.plot)+5)

dev.off()
##########################################################################################################################################

saveWorkbook(wb, file = paste0(panel.directory,"Source_data_Fig.3.xlsx"), overwrite=T)
saveWorkbook(wb.s2, file = paste0(output.directory, "Figure2/Source_data_Fig.S2.xlsx"), overwrite=T)
saveWorkbook(wb.s3, file = paste0(output.directory, "Figure2/Source_data_Fig.S3.xlsx"), overwrite=T)


