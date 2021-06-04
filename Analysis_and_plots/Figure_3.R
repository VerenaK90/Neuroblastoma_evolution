## Reproduce Fig. 3
##############################################################################################################################################
## load settings and libraries

load(paste0(rdata.directory, "MRCA_timing.RData"))

## source data:
wb <- createWorkbook()

## store figure panels
panel.directory <- paste0(output.directory, "Figure3/")

if(!dir.exists(panel.directory)){
  dir.create(panel.directory)
}



##############################################################################################################################################
## Figure 3a-d: validation cohort: density distributions for MRCA and ECA 

source("./Nextcloud/NB_manuscript/Revised_version/Plots_and_scripts/Custom_scripts/Oncoprint_additional_samples.R")

max.mutation.time.primary <- max(mutation.time.mrca[rownames(sample.information.30x[sample.information.30x$Location %in% c("Primary", "Metastasis"),])],
                                 na.rm=T)
  
sample.information.30x$Telomere.maintenance.mechanism <- factor(sample.information.30x$Telomere.maintenance.mechanism,
                                                        levels=c("MNA", "TERT", "ALT", "Multiple", "None"))

##### Primary tumors, with and without treatment, metastases
## Plot the cases with ECA separately per ploidy; merge di- and the tetraploid case.

p1 <- list()
p1.subset.list <- list()

for(i in c("triploid" ,"di/tetraploid")){
  
  if(i=="di/tetraploid"){
    i <- c(2,4)
  }else{
    i <- 3
  }
  print(i)
  
  subset=sample.information.30x[sample.information.30x$Rounded.ploidy %in% i & sample.information.30x$Location %in% c("Primary", "Metastasis") &
                           sample.information.30x$ECA.exists=="T",,drop=F]
  
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
  
  
  
  if(i%in%c(2,4)){
    panel="b"
    addWorksheet(wb, panel)
  }else if(i==3){
    panel="a"
    addWorksheet(wb, panel)
  }
  
  
  writeData(wb, panel, to.plot[,c("ECA", "ECA.lower", "ECA.upper", "MRCA", "MRCA.lower", "MRCA.upper")], )
  
  
  
  p1[[length(p1)+1]] <- ggplot(data = to.plot[order(to.plot$MRCA),],
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
                    fill=unname(manual.colors["Early"]), alpha=0.5, col=NA) +
    scale_x_continuous(name = "Mutations/Mb",
                       limits=c(0, max.mutation.time.primary/3.3/10^3))+
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(name = "Fraction of tumors") +
    ggtitle(paste("Ploidy =",i, " # Cases = ", nrow(subset)))
  
  ## TMM, Clinical subtype & CNVs
  to.plot. <- subset[,c("MRCAtime", "Telomere.maintenance.mechanism", "ManualScore")]
  to.plot.$xmin <-c(0, (to.plot.$MRCAtime/3.3/10^3)[-length(to.plot.$MRCAtime)])
  to.plot.$xmax <- to.plot.$MRCAtime/3.3/10^3
  to.plot.$ymin <- 0
  to.plot.$ymax <- 0.23
  
  to.plot.[,c("17q", "17", "1p", "1q", "1", "7q", "7", "2p", "2", "11q")] <- matrix.early.late.30x[c("17q", "17", "1p", "1q", "1", "7q", "7", "2p", "2", "11q"),
                                                                                               rownames(to.plot.)]
  
  ## distinguish cases that are not dateable from cases that are statistically insigniicant from ECA
  for(chr.change in c("17q", "17", "1p", "1q", "1", "7q", "7", "2p", "2", "11q")){
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
  
  p1.subset.list[[length(p1.subset.list) + 1]] <- ggplot(data=to.plot.,
                                                         aes(xmin=xmin, xmax=xmax,
                                                             ymin=ymin, ymax=ymax,
                                                             fill=as.character(Telomere.maintenance.mechanism))) +
    geom_rect() +
    geom_rect(data=to.plot., aes(xmin= xmin, xmax=xmax, ymin=ymin+0.27, ymax=ymin+0.48, fill=ManualScore), inherit.aes = F) +
    scale_fill_manual(values=c(telomere.colors, clinical.risk.colors, manual.colors, "n.d."="grey")) +
    scale_x_continuous(limits=c(0, max.mutation.time.primary/3.3/10^3)) +
    theme(legend.position = "bottom") +
    scale_y_continuous(breaks=seq(0.125, 2.875, 0.25),
                       labels=c("Stage", "TMM", rev(c("17q", "17", "1p", "1q", "1", "7q", "7", "2p", "2", "11q"))))
  
  for(chr.change in rev(c("`17q`", "`17`", "`1p`", "`1q`", "`1`", "`7q`", "`7`", "`2p`", "`2`", "`11q`"))){
    
    index <- which(rev(c("`17q`", "`17`", "`1p`", "`1q`", "`1`", "`7q`", "`7`", "`2p`", "`2`", "`11q`"))==chr.change)
    to.plot.$ymin <- 0.5 + 0.25*(index - 1)+0.2
    to.plot.$ymax <- 0.5 + 0.25*index-0.2
    
    p1.subset.list[[length(p1.subset.list)]] <- p1.subset.list[[length(p1.subset.list)]] +
      geom_rect(data=to.plot., aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill=chr.change), inherit.aes=T)
  }
  
  
  
}


## Plot the cases without ECA separately per ploidy; merge di- and tetraploids

p2 <- list()
p2.subset.list <- list()

for(i in c("triploid", "di/tetraploid")){
  
  if(i=="di/tetraploid"){
    i <- c(2,4)
  }else{
    i <- 3
  }
  
  
  
  subset=sample.information.30x[sample.information.30x$Rounded.ploidy %in% i & sample.information.30x$Location %in% c("Primary", "Metastasis") &
                           sample.information.30x$ECA.exists=="F",,drop=F]
  
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
  
  
  ## export the source data
  if(i==2){
    
    
    panel="d"
    
    addWorksheet(wb, panel)
    writeData(wb, panel, to.plot[,c("ECA", "ECA.lower", "ECA.upper", "MRCA", "MRCA.lower", "MRCA.upper")])
    
  }else{
    panel="c"
    
    addWorksheet(wb, panel)
    writeData(wb, panel, to.plot[,c("ECA", "ECA.lower", "ECA.upper", "MRCA", "MRCA.lower", "MRCA.upper")])
    
  }
  
  
  p2[[length(p2)+1]] <- ggplot(data = to.plot[order(to.plot$MRCA),],
                               aes(x=MRCA, y=seq(1/length(MRCA),1,length.out = length(MRCA)),
                                   ymin =  sapply(sort(MRCA), function(x){
                                     sum(MRCA.upper <= x)
                                   })/length(MRCA),
                                   ymax= sapply(sort(MRCA), function(x){
                                     sum(MRCA.lower <= x)
                                   })/length(MRCA)
                               )) +
    stat_ecdf(col=unname(manual.colors["Late"])) +
    geom_stepribbon(fill=unname(manual.colors["Late"]), alpha=0.5, col=NA) +
    scale_x_continuous(name = "Mutations/Mb", limits=c(0, max.mutation.time.primary/3.3/10^3))+
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(name = "Fraction of tumors") +
    ggtitle(paste("Ploidy =",i, ", # Cases = ", nrow(subset)))
  
  
  ## TMM, Clinical subtype & CNVs
  to.plot. <- subset[,c("MRCAtime", "Telomere.maintenance.mechanism", "ManualScore")]
  to.plot.$xmin <-c(0, (to.plot.$MRCAtime/3.3/10^3)[-length(to.plot.$MRCAtime)])
  to.plot.$xmax <- to.plot.$MRCAtime/3.3/10^3
  to.plot.$ymin <- 0.02
  to.plot.$ymax <- 0.23
  
  to.plot.[,c("17q", "17", "1p", "1q", "1", "7q", "7", "2p", "2", "11q")] <- t(matrix.early.late.30x[c("17q", "17", "1p", "1q", "1", "7q", "7", "2p", "2", "11q"),
                                                                                                 rownames(to.plot.)])
  
  ## distinguish cases that are not dateable from cases that are statistically insigniicant from ECA
  for(chr.change in c("17q", "17", "1p", "1q", "1", "7q", "7", "2p", "2", "11q")){
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
  
  p1.subset.list[[length(p1.subset.list) + 1]] <- ggplot(data=to.plot.,
                                                         aes(xmin=xmin, xmax=xmax,
                                                             ymin=ymin, ymax=ymax,
                                                             fill=as.character(Telomere.maintenance.mechanism))) +
    geom_rect() +
    geom_rect(data=to.plot., aes(xmin= xmin, xmax=xmax, ymin=ymin+0.27, ymax=ymin+0.48, fill=ManualScore), inherit.aes = F) +
    scale_fill_manual(values=c(telomere.colors, clinical.risk.colors, manual.colors, "n.d." = "grey")) +
    scale_x_continuous(limits=c(0, max.mutation.time.primary/3.3/10^3)) +
    theme(legend.position = "bottom") +
    scale_y_continuous(breaks=seq(0.125, 2.875, 0.25),
                       labels=c("Stage", "TMM", rev(c("17q", "17", "1p", "1q", "1", "7q", "7", "2p", "2", "11q"))))
  
  for(chr.change in rev(c("`17q`", "`17`", "`1p`", "`1q`", "`1`", "`7q`", "`7`", "`2p`", "`2`", "`11q`"))){
    
    index <- which(rev(c("`17q`", "`17`", "`1p`", "`1q`", "`1`", "`7q`", "`7`", "`2p`", "`2`", "`11q`"))==chr.change)
    to.plot.$ymin <- 0.5 + 0.25*(index - 1) + 0.02
    to.plot.$ymax <- 0.5 + 0.25*index - 0.02
    
    p1.subset.list[[length(p1.subset.list)]] <- p1.subset.list[[length(p1.subset.list)]] +
      geom_rect(data=to.plot., aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill=chr.change), inherit.aes=T)
  }
  
}


pdf(paste0(panel.directory, "Figure_3a_b_c_d.pdf"), width=9, height=9, useDingbats = F)

figure <- ggarrange(plotlist=c(p1,p2), nrow=3, ncol=3)
annotate_figure(figure, top="Primary tumor / Metastasis")

figure <- ggarrange(plotlist=c(p1.subset.list,p2.subset.list), nrow=3, ncol=3)
annotate_figure(figure, top="Primary tumor / Metastasis")

dev.off()


##############################################################################################################################################
saveWorkbook(wb, file = paste0(panel.directory,"Source_data_Fig.3.xlsx"), overwrite=T)


##############################################################################################################################################
## Some statistical tests:

## mean and sd mutation density at triploidization
mutation.time.eca[is.na(mutation.time.eca)] <- mutation.time.mrca[is.na(mutation.time.eca)]
## all 30x tumors are primaries
mean(mutation.time.eca[c(intersect(primary.tumors.80x, triploid.tumors.80x), triploid.tumors.30x)]/3.3/10^3)
sd(mutation.time.eca[c(intersect(primary.tumors.80x, triploid.tumors.80x), triploid.tumors.30x)]/3.3/10^3)


## compare mean and sd mutation density at individual chromosomal gains/genome doubling between datasets 
## all 30x tumors are primaries
mean(mutation.time.eca[c(intersect(primary.tumors.80x, c(diploid.tumors.80x, tetraploid.tumors.80x)),
                         diploid.tumors.30x, tetraploid.tumors.30x)]/3.3/10^3)
sd(mutation.time.eca[c(intersect(primary.tumors.80x, c(diploid.tumors.80x, tetraploid.tumors.80x)),
                       diploid.tumors.30x, tetraploid.tumors.30x)]/3.3/10^3)


## mutation density of MRCAs with single hits
mean(mutation.time.mrca[mutation.time.mrca==mutation.time.eca & 
                          names(mutation.time.mrca) %in% c(intersect(primary.tumors.80x, c(diploid.tumors.80x, tetraploid.tumors.80x)),
                                                           diploid.tumors.30x, tetraploid.tumors.30x)])/3.3/10^3
sd(mutation.time.mrca[mutation.time.mrca==mutation.time.eca & 
                        names(mutation.time.mrca) %in% c(intersect(primary.tumors.80x, c(diploid.tumors.80x, tetraploid.tumors.80x)),
                                                         diploid.tumors.30x, tetraploid.tumors.30x)])/3.3/10^3

mean(mutation.time.mrca[mutation.time.mrca==mutation.time.eca & 
                          names(mutation.time.mrca) %in% c(intersect(primary.tumors.80x, triploid.tumors.80x), triploid.tumors.30x)])/3.3/10^3
sd(mutation.time.mrca[mutation.time.mrca==mutation.time.eca & 
                          names(mutation.time.mrca) %in% c(intersect(primary.tumors.80x, triploid.tumors.80x), triploid.tumors.30x)])/3.3/10^3
