## Reproduce Fig. 3
##############################################################################################################################################
## load settings and libraries

load(paste0(rdata.directory, "MRCA_timing.RData"))

## source data:
wb <- createWorkbook()
## add to Supplementary Fig. 2
wb.s <- loadWorkbook(paste0(paste0(output.directory, "Figure2/"), "Source_data_Fig.S2.xlsx"))

## store figure panels
panel.directory <- paste0(output.directory, "Figure3/")

if(!dir.exists(panel.directory)){
  dir.create(panel.directory)
}


##############################################################################################################################################
## Figure 3a: Mutation density at MRCAs in the discovery set
## color by cutpoint
cutpoint <- 0.05

pdf(paste0(panel.directory, "Figure_3a.pdf"), width=3.5, height=3.5, useDingbats = F)

to.plot <- data.frame(MRCA=mutation.time.mrca[primary.tumors.80x])
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
## Figure 3b,c: density distributions for MRCA and ECA - discovery set
source(paste0(custom.script.directory, "Driver_mutations.R"))
source(paste0(custom.script.directory, "Oncoprint.R"))

## if ECA was not uniquely identified, take the earliest time point
mutation.time.eca[names(earliest.mutation.time)] <- earliest.mutation.time
mutation.time.eca.lower[names(earliest.mutation.time)] <- earliest.mutation.time.lower
mutation.time.eca.upper[names(earliest.mutation.time)] <- earliest.mutation.time.upper
####

colnames(sample.information.80x)[which(colnames(sample.information.80x)=="ECA")] <- "ECA.exists"

max.mutation.time.primary <- max(mutation.time.mrca[rownames(sample.information.80x[sample.information.80x$Location %in% c("Primary", "Relapse"),])],
                                 na.rm=T)

sample.information.80x$Telomere.maintenance.mechanism <- factor(sample.information.80x$Telomere.maintenance.mechanism,
                                                                levels=c("MNA", "TERT", "ALT", "Multiple", "None"))

sample.information.80x$ECA.exists <- as.logical(sample.information.80x$ECA.exists)

##### Primary tumors, with and without treatment, metastases
## Plot early and late MRCA separately; further, stratify late cases w.r.t. number of events

cutpoint <- 0.05

plate <- list()
plate.subset.list <- list()

for(ECA.exists in c(T, F)){
  
  subset=sample.information.80x[ sample.information.80x$Location %in% c("Primary", "Metastasis") &
                                   sample.information.80x$ECA.exists==ECA.exists &
                                   mutation.time.mrca[rownames(sample.information.80x)]/3.3/10^3 >=cutpoint,,drop=F]
  
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
    panel="c_ECA"
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
    panel="c_ECA_sample_info"
  }else{
    panel="c_no_ECA_sample_info"
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
  
  subset=sample.information.80x[ sample.information.80x$Location %in% c("Primary", "Metastasis") &
                                   mutation.time.mrca[rownames(sample.information.80x)]/3.3/10^3 < cutpoint,,drop=F]
  
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
  
  panel="b_sample_info"
  
  addWorksheet(wb, panel)
  to.plot.$Stage <- replace(to.plot.$ManualScore, to.plot.$ManualScore=="LR", "1, 2, 4S")
  to.plot.$Stage <- replace(to.plot.$Stage, to.plot.$Stage=="IR", "3")
  to.plot.$Stage <- replace(to.plot.$Stage, to.plot.$Stage=="HR", "4")
  writeData(wb, panel, to.plot.[,c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2",
                                   "Telomere.maintenance.mechanism", "Stage", "Ploidy")]) 
  
}


pdf(paste0(panel.directory, "Figure_3_b_c.pdf"), width=9, height=9, useDingbats = F)

figure <- ggarrange(plotlist=c(pearly.unsplit, plate), nrow=3, ncol=3)
annotate_figure(figure, top="Primary tumor / Metastasis")

figure <- ggarrange(plotlist=c(pearly.unsplit.subset.list, plate.subset.list), nrow=3, ncol=3)
annotate_figure(figure, top="Primary tumor / Metastasis")

dev.off()


##########################################################################################################################################
## Figure S2h: Compare MRCA between cases with and without ECA among late-MRCA tumors

pdf(paste0(panel.directory,"Figure_S2h.pdf"), useDingbats = F)


to.plot <- data.frame(MRCA=mutation.time.mrca[primary.tumors.80x]/3.3/10^3,
                      ECA.exists=ifelse(is.na(mutation.time.eca[primary.tumors.80x]) | mutation.time.eca[primary.tumors.80x]!=mutation.time.mrca[primary.tumors.80x], T, F),
                      MRCA.time=ifelse(mutation.time.mrca[primary.tumors.80x]/3.3/10^3<cutpoint, "early", "late"))

addWorksheet(wb.s, "2h")
writeData(wb.s, "2h", to.plot)

ggplot(to.plot[to.plot$MRCA.time=="late" ,], 
       aes(x=ECA.exists, y=MRCA, col=ECA.exists)) + geom_boxplot() + geom_beeswarm() +
  scale_color_manual(values=c("TRUE" = unname(manual.colors["Early"]), "FALSE" = unname(manual.colors["Late"]))) +
  geom_boxplot(data=to.plot[to.plot$MRCA.time=="early",], col="dodgerblue", aes(x="Early MRCA"))+
  geom_beeswarm(data=to.plot[to.plot$MRCA.time=="early",], col="dodgerblue", aes(x="Early MRCA"))


wilcox.test(to.plot[to.plot$MRCA.time=="late" & to.plot$ECA.exists==T,]$MRCA, 
            to.plot[to.plot$MRCA.time=="late" & to.plot$ECA.exists==F,]$MRCA)

## p = 6.4 e-05


dev.off()

##########################################################################################################################################
#### Figure S2i: # of events for triploidization

mutation.time.eca[is.na(mutation.time.eca)] <- mutation.time.mrca[is.na(mutation.time.eca)]

## For each ploidy: can we assign most polyploidization events to a single event? 
numbers.of.3n.chromosomes.at.single.event <- c()
numbers.of.3n.chromosomes <- c()
not.conforming.chromosomes <- list()

for(i in triploid.tumors.80x){
  
  ## take only dominating copy number per chromosome
  
  dominating.copy.number.per.chromosome <- paste("chr", 1:22, lapply(segment.length.matrix, function(x){
    rownames(x[-c(4,7),])[which.max(x[-c(4,7),i])]
  }), sep="_")
  
  if(mutation.time.eca[i]==mutation.time.mrca[i]){
    
    ## if these have a higher mutation count than the mrca estimate count them to mrca
    tmp <- c(gains.at.mrca[[i]], gains.not.maping.to.eca.or.mrca[[i]][mutation.time.most.likely[[i]][gains.not.maping.to.eca.or.mrca[[i]]] >=
                                                                        mutation.time.mrca[i]])
    tmp <- sapply(intersect(tmp, dominating.copy.number.per.chromosome), function(x){strsplit(x, split="_")[[1]][3]})
    
    tmp.2 <- gains.not.maping.to.eca.or.mrca[[i]][mutation.time.most.likely[[i]][gains.not.maping.to.eca.or.mrca[[i]]] <
                                                    mutation.time.mrca[i]]   
    tmp.2 <- intersect(tmp.2, dominating.copy.number.per.chromosome)
    tmp.2 <- sapply(intersect(tmp.2, dominating.copy.number.per.chromosome), function(x){strsplit(x, split="_")[[1]][3]})
    
    
    numbers.of.3n.chromosomes.at.single.event[i] <- sum(tmp %in% c("trisomic", "tetrasomic"))
    numbers.of.3n.chromosomes[i] <- (sum(tmp.2 %in% c("trisomic", "tetrasomic"))+sum(tmp %in% c("trisomic", "tetrasomic")))
    
    tmp.2 <- sapply(names(tmp.2[tmp.2 %in% c("trisomic", "tetrasomic")]), function(x){strsplit(x, split="_")[[1]][2]})
    not.conforming.chromosomes[[i]] <- tmp.2
    
  }else{
    tmp <- c(gains.uniquely.mapped.to.eca[[i]], gains.at.mrca.conforming.eca[[i]])
    tmp <- sapply(intersect(tmp, dominating.copy.number.per.chromosome), function(x){strsplit(x, split="_")[[1]][3]})
    
    tmp.2 <- c(setdiff(c(gains.at.mrca[[i]], gains.not.maping.to.eca.or.mrca[[i]][mutation.time.most.likely[[i]][gains.not.maping.to.eca.or.mrca[[i]]] >=
                                                                                    mutation.time.mrca[i]]), gains.at.mrca.conforming.eca[[i]]), 
               gains.not.maping.to.eca.or.mrca[[i]][mutation.time.most.likely[[i]][gains.not.maping.to.eca.or.mrca[[i]]] < mutation.time.mrca[i]])
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

for(i in tetraploid.tumors.80x){
  
  ## take only dominating copy number per chromosome
  
  dominating.copy.number.per.chromosome <- paste("chr", 1:22, lapply(segment.length.matrix, function(x){
    rownames(x[-c(4,7),])[which.max(x[-c(4,7),i])]
  }), sep="_")
  
  if(mutation.time.eca[i]==mutation.time.mrca[i]){
    
    ## if these have a higher mutation count than the mrca estimate count them to mrca
    tmp <- c(gains.at.mrca[[i]], gains.not.maping.to.eca.or.mrca[[i]][mutation.time.most.likely[[i]][gains.not.maping.to.eca.or.mrca[[i]]] >=
                                                                        mutation.time.mrca[i]])
    tmp <- sapply(intersect(tmp, dominating.copy.number.per.chromosome), function(x){strsplit(x, split="_")[[1]][3]})
    
    ## if these have a higher mutation count than the mrca estimate count them to mrca
    tmp.2 <- gains.not.maping.to.eca.or.mrca[[i]][mutation.time.most.likely[[i]][gains.not.maping.to.eca.or.mrca[[i]]] <
                                                    mutation.time.mrca[i]]
    
    tmp.2 <- intersect(tmp.2, dominating.copy.number.per.chromosome)
    tmp.2 <- sapply(intersect(tmp.2, dominating.copy.number.per.chromosome), function(x){strsplit(x, split="_")[[1]][3]})
    
    
    numbers.of.4n.chromosomes.at.single.event[i] <- sum(tmp %in% c("trisomic", "tetrasomic"))
    numbers.of.4n.chromosomes[i] <- (sum(tmp.2 %in% c("trisomic", "tetrasomic"))+sum(tmp %in% c("trisomic", "tetrasomic")))
    
    tmp.2 <- sapply(names(tmp.2[tmp.2 %in% c("trisomic", "tetrasomic")]), function(x){strsplit(x, split="_")[[1]][2]})
    not.conforming.chromosomes[[i]] <- tmp.2
    
  }else{
    tmp <- c(gains.uniquely.mapped.to.eca[[i]], gains.at.mrca.conforming.eca[[i]])
    tmp <- sapply(intersect(tmp, dominating.copy.number.per.chromosome), function(x){strsplit(x, split="_")[[1]][3]})
    
    tmp.2 <- c(setdiff(c(gains.at.mrca[[i]], gains.not.maping.to.eca.or.mrca[[i]][mutation.time.most.likely[[i]][gains.not.maping.to.eca.or.mrca[[i]]] >=
                                                                                    mutation.time.mrca[i]]), gains.at.mrca.conforming.eca[[i]]), 
               gains.not.maping.to.eca.or.mrca[[i]][mutation.time.most.likely[[i]][gains.not.maping.to.eca.or.mrca[[i]]] < mutation.time.mrca[i]])
    
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

for(i in diploid.tumors.80x){
  
  ## take only dominating copy number per chromosome
  
  dominating.copy.number.per.chromosome <- paste("chr", 1:22, lapply(segment.length.matrix, function(x){
    rownames(x[-c(4,7),])[which.max(x[-c(4,7),i])]
  }), sep="_")
  
  ## among them select the disomic chromosomes
  
  disomic.chromosomes <- dominating.copy.number.per.chromosome[sapply(dominating.copy.number.per.chromosome, function(x){
    x <- strsplit(x, split="_")[[1]]
    if(x[3]=="disomic"){
      return(T)
    }else{
      return(F)
    }
  })]
  
  tmp <- setdiff(disomic.chromosomes, monosomic.states.not.matching.mrca[[i]])
  if(length(monosomic.states.not.matching.mrca[[i]])==0){
    tmp.2 <- monosomic.states.not.matching.mrca[[i]]
  }else{
    tmp.2 <- monosomic.states.not.matching.mrca[[i]][sapply(monosomic.states.not.matching.mrca[[i]], function(x){
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

pdf(paste0(panel.directory,"Figure_S2i.pdf"), width=2, height=2.2, useDingbats = F)

to.plot <- rbind(data.frame(Chromosomes=numbers.of.2n.chromosomes.at.mrca/numbers.of.2n.chromosomes,
                            Tumor=diploid.tumors.80x,
                            Telomere.type=telomere.classification.80x[diploid.tumors.80x],
                            Time = sample.information.80x[diploid.tumors.80x,]$Location,
                            Ploidy=2),
                 data.frame(Chromosomes=numbers.of.3n.chromosomes.at.single.event/numbers.of.3n.chromosomes,
                            Tumor=triploid.tumors.80x,
                            Telomere.type=telomere.classification.80x[triploid.tumors.80x],
                            Time = sample.information.80x[triploid.tumors.80x,]$Location,
                            Ploidy=3),
                 data.frame(Chromosomes=numbers.of.4n.chromosomes.at.single.event/numbers.of.4n.chromosomes,
                            Tumor=tetraploid.tumors.80x,
                            Telomere.type=telomere.classification.80x[tetraploid.tumors.80x],
                            Time = sample.information.80x[tetraploid.tumors.80x,]$Location,
                            Ploidy=4))


panel="i"
addWorksheet(wb.s, panel)
writeData(wb.s, panel, to.plot)


ggplot(to.plot, aes(x=Ploidy, y=Chromosomes,group=Ploidy)) + geom_boxplot(width=0.5) + geom_beeswarm() + scale_y_continuous(limits=c(0,1.01))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


dev.off()


##########################################################################################################################################
## Fig. 3d: compare mutation density at ECA of late-MRCA group with MRCA of early-MRCAgroup

pdf(paste0(panel.directory,"Figure_3d.pdf"), useDingbats = F, width=2, height=2)

to.plot <- data.frame(MRCA=mutation.time.mrca[primary.tumors.80x]/3.3/10^3,
                      ECA=mutation.time.eca[primary.tumors.80x]/3.3/10^3,
                      ECA.exists=ifelse(is.na(mutation.time.eca[primary.tumors.80x]) | mutation.time.eca[primary.tumors.80x]!=mutation.time.mrca[primary.tumors.80x], T, F),
                      MRCA.time=ifelse(mutation.time.mrca[primary.tumors.80x]/3.3/10^3<cutpoint, "early", "late"))

to.plot <- to.plot[to.plot$ECA.exists==T | to.plot$MRCA.time=="early",]

addWorksheet(wb, "3d")
writeData(wb, "3d", to.plot)

ggplot(to.plot[to.plot$MRCA.time=="late" ,], 
       aes(x=ECA.exists, y=ECA)) + geom_boxplot(col="dodgerblue4") + geom_beeswarm(col="dodgerblue4") +
  geom_boxplot(data=to.plot[to.plot$MRCA.time=="early",], col="dodgerblue", aes(x="Early MRCA", y=MRCA))+
  geom_beeswarm(data=to.plot[to.plot$MRCA.time=="early",], col="dodgerblue", aes(x="Early MRCA", y=MRCA)) +
  scale_y_continuous(name="#SSNVs/Mb")

wilcox.test(to.plot[to.plot$MRCA.time=="late",]$ECA, 
            to.plot[to.plot$MRCA.time=="early",]$MRCA)

## p=0.23

dev.off()


##########################################################################################################################################

saveWorkbook(wb, file = paste0(panel.directory,"Source_data_Fig.3.xlsx"), overwrite=T)
saveWorkbook(wb.s, file = paste0(output.directory, "Figure2/Source_data_Fig.S2.xlsx"), overwrite=T)
