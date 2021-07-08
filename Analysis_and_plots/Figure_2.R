## Reproduce Fig. 2
##############################################################################################################################################
## load settings and libraries

## either re-run
source(paste0(custom.script.directory, "Mutation_density_quantification.R"))
## or load
#load(paste0(rdata.directory, "Clonal_mutations_different_ploidies.RData"))
#load(paste0(rdata.directory, "MRCA_timing.RData"))
#load(paste0(rdata.directory, "Vafs_all_tumors.RData"))

## get mutation rate per day from later analysis (Figure 5; code this differently/in order)

## source data:
wb <- createWorkbook()
wb.s <- createWorkbook()

## store figure panels
panel.directory <- paste0(output.directory, "Figure2/")

if(!dir.exists(panel.directory)){
  dir.create(panel.directory)
}

##############################################################################################################################################
## Figure 2b: NBE15 

id <- "NBE15"

to.plot <- data.frame(VAF = c(vafs.all.tumors[[id]][[2]][,2]/rowSums(vafs.all.tumors[[id]][[2]]),
                      vafs.all.tumors[[id]][[3]][,2]/rowSums(vafs.all.tumors[[id]][[3]]),
                      vafs.all.tumors[[id]][[4]][,2]/rowSums(vafs.all.tumors[[id]][[4]])),
                      Ploidy = c(rep(2, nrow(vafs.all.tumors[[id]][[2]])),
                                 rep(3, nrow(vafs.all.tumors[[id]][[3]])),
                                 rep(4, nrow(vafs.all.tumors[[id]][[4]]))))
addWorksheet(wb, "b")
writeData(wb, "b", to.plot)

pdf(paste0(panel.directory, "Figure_2b.pdf"))

ggplot(to.plot[to.plot$Ploidy==2,], aes(x=VAF)) + geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(0,1), name="VAF") +
  ggtitle(label ="Disomic chromosomes")
ggplot(to.plot[to.plot$Ploidy==3,], aes(x=VAF)) + geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(0,1), name="VAF") +
  ggtitle(label="Trisomic chromosomes")
ggplot(to.plot[to.plot$Ploidy==4,], aes(x=VAF)) + geom_histogram(binwidth=0.01) + scale_x_continuous(limits=c(0,1), name="VAF") +
  ggtitle(label="Tetrasomic chromosomes")

dev.off()


##############################################################################################################################################
## Figure 2d-f and Fig. S2c,d,e: density plots per genomic segment for NBE15, NBE14, NBE44, NBE99


## manually set the range 
max.mutation.time <- 0.7

## workbook for source data

for(i in c("NBE15", "NBE14", "NBE44", "NBE99")){
  if(length( mutation.time.most.likely[[i]])==0){next}
  
  binwidth = (max(mutation.time.most.likely[[i]]/3.3/10^3) - min(mutation.time.most.likely[[i]]/3.3/10^3))/20
  
  to.plot <- data.frame(x= mutation.time.most.likely[[i]]/3.3/10^3,
                        xupper= mutation.time.upper[[i]]/3.3/10^3,
                        xlower=mutation.time.lower[[i]]/3.3/10^3,
                        Segment = names(mutation.time.most.likely[[i]]),
                        Timing = sapply(names(mutation.time.most.likely[[i]]), function(x){
                          x <- strsplit(x, split="_")[[1]]
                          if(is.na(x[4])){
                            return("Post-CNV")
                          }
                          if(x[4]=="I"){
                            return("Post-CNV")
                          }else{
                            return("Pre-CNV")
                          }
                        }))
  
  to.plot$Segment <- factor(to.plot$Segment, levels=unique(to.plot$Segment))
  ## for the three examples shown in the manuscript, export the source data
  if(i %in% c("NBE15", "NBE14", "NBE44", "NBE99")){
    
    if(i=="NBE15"){
      panel="d_e_f"
      addWorksheet(wb, panel)
      writeData(wb, panel, to.plot)
    }else if(i=="NBE14"){
      panel="d"
      addWorksheet(wb.s, panel)
      writeData(wb.s, panel, to.plot)
    }else if(i=="NBE99"){
      panel="c"
      addWorksheet(wb.s, panel)
      writeData(wb.s, panel, to.plot)
    }else{
      panel="e"
      addWorksheet(wb.s, panel)
      writeData(wb.s, panel, to.plot)
    }
    

    
  }
  
 
  if(i == "NBE15"){
    pdf(paste0(panel.directory, "Figure_2d_e_f.pdf"))
    
    p <- ggplot(to.plot[to.plot$Timing=="Post-CNV",], aes(x=x)) + geom_histogram( binwidth = binwidth, fill=manual.colors["Late"])+ 
      scale_y_continuous(name="# Genomic segments") + 
       scale_x_continuous(name="Mutations per Mb", limits = c(-0.05*(max(to.plot$x)), max.mutation.time)) +
      ggtitle(i)
    
    print(p)
    
    p <- ggplot(to.plot[to.plot$Timing=="Pre-CNV",], aes(x=x, fill=Segment)) + geom_histogram( binwidth = binwidth)+ 
      # scale_fill_manual(values=brewer.pal(9, "PRGn")[c(9)])+
      scale_fill_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot[to.plot$Timing=="Pre-CNV",])))+
      geom_density(data=to.plot[to.plot$Timing=="Post-CNV",], linetype=2, aes(x=x), inherit.aes = F) +
      scale_y_continuous(name="# Genomic segments") + 
        scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$x)), max(to.plot$x)*1.1)) +
      ggtitle(i)
    
    print(p)
    
    to.plot. <- to.plot[to.plot$Timing=="Pre-CNV",]
    to.plot.$Segment <- factor(to.plot.$Segment, levels=unique(to.plot.$Segment))
    
    p <- ggplot(to.plot., aes(x=x, xmin=xlower, xmax=xupper, y=as.numeric(Segment)/nrow(to.plot.), col=Segment)) +
      scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot.))) + 
      geom_point() + geom_errorbarh(height=0)+ 
      geom_vline(data=data.frame(x=mutation.time.mrca[i]), aes(xintercept=x/3.3/10^3, y=1), inherit.aes = F, col=manual.colors["Late"],
                 linetype=2) + 
      geom_ribbon(data=data.frame(xmin=rep(mutation.time.mrca.lower[i],2),
                                  xmax=rep(mutation.time.mrca.upper[i],2),
                                  y=c(0, nrow(to.plot[to.plot$Timing=="Pre-CNV",]))), 
                  aes(xmin=xmin/3.3/10^3, xmax=xmax/3.3/10^3,y=y),
                  inherit.aes = F, fill=manual.colors["Late"], col=NA, alpha=0.5) +
      geom_vline(data=data.frame(x=mutation.time.eca[i]), aes(xintercept=x/3.3/10^3, y=1), inherit.aes = F, col=manual.colors["Early"],
                 linetype=2) + 
      geom_ribbon(data=data.frame(xmin=rep(mutation.time.eca.lower[i], 2),
                                  xmax=rep(mutation.time.eca.upper[i],2),
                                  y=c(0, nrow(to.plot[to.plot$Timing=="Pre-CNV",]))), 
                  aes( xmin=xmin/3.3/10^3, xmax=xmax/3.3/10^3, y=y),
                  height=0, inherit.aes = F, fill=manual.colors["Early"], col=NA, alpha=0.5) +
      scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$x)), max(to.plot$xupper)*1.1)) +
      scale_y_continuous(limits=c(0, nrow(to.plot)))
    
    print(p)

    dev.off()
    
  }
  
  if(i == "NBE14"){
    pdf(paste0(panel.directory, "Figure_S2d.pdf"))
    
    p <- ggplot(to.plot[to.plot$Timing=="Post-CNV",], aes(x=x)) + geom_histogram( binwidth = binwidth, fill=manual.colors["Late"])+ 
      scale_y_continuous(name="# Genomic segments") + 
        scale_x_continuous(name="Mutations per Mb", limits = c(-0.05*(max(to.plot$x)), max.mutation.time)) +
      ggtitle(i)
    
    print(p)
    
    p <- ggplot(to.plot[to.plot$Timing=="Pre-CNV",], aes(x=x, fill=Segment)) + geom_histogram( binwidth = binwidth)+ 
      scale_fill_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot[to.plot$Timing=="Pre-CNV",])))+
      geom_density(data=to.plot[to.plot$Timing=="Post-CNV",], linetype=2, aes(x=x), inherit.aes = F) +
      scale_y_continuous(name="# Genomic segments") + 
        scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$x)), max(to.plot$x)*1.1)) +
      ggtitle(i)
    
    print(p)
    
    to.plot. <- to.plot[to.plot$Timing=="Pre-CNV",]
    to.plot.$Segment <- factor(to.plot.$Segment, levels=unique(to.plot.$Segment))
    
    p <- ggplot(to.plot., aes(x=x, xmin=xlower, xmax=xupper, y=as.numeric(Segment)/nrow(to.plot.), col=Segment)) +
      scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot.))) + 
      geom_point() + geom_errorbarh(height=0)+ 
      geom_vline(data=data.frame(x=mutation.time.mrca[i]), aes(xintercept=x/3.3/10^3, y=1), inherit.aes = F, col=manual.colors["Late"],
                 linetype=2) + 
      geom_ribbon(data=data.frame(xmin=rep(mutation.time.mrca.lower[i],2),
                                  xmax=rep(mutation.time.mrca.upper[i],2),
                                  y=c(0, nrow(to.plot[to.plot$Timing=="Pre-CNV",]))), 
                  aes(xmin=xmin/3.3/10^3, xmax=xmax/3.3/10^3,y=y),
                  inherit.aes = F, fill=manual.colors["Late"], col=NA, alpha=0.5) +
      geom_vline(data=data.frame(x=mutation.time.eca[i]), aes(xintercept=x/3.3/10^3, y=1), inherit.aes = F, col=manual.colors["Early"],
                 linetype=2) + 
      geom_ribbon(data=data.frame(xmin=rep(mutation.time.eca.lower[i], 2),
                                  xmax=rep(mutation.time.eca.upper[i],2),
                                  y=c(0, nrow(to.plot[to.plot$Timing=="Pre-CNV",]))), 
                  aes( xmin=xmin/3.3/10^3, xmax=xmax/3.3/10^3, y=y),
                  height=0, inherit.aes = F, fill=manual.colors["Early"], col=NA, alpha=0.5) +
      scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$x)), max(to.plot$xupper)*1.1)) +
      scale_y_continuous(limits=c(0, nrow(to.plot)))
    
    print(p)

    dev.off()
  }
  
  if(i == "NBE44"){
    pdf(paste0(panel.directory, "Figure_S2e.pdf"))
    
    p <- ggplot(to.plot[to.plot$Timing=="Post-CNV",], aes(x=x)) + geom_histogram( binwidth = binwidth, fill=manual.colors["Late"])+ 
      scale_y_continuous(name="# Genomic segments") + 
       scale_x_continuous(name="Mutations per Mb", limits = c(-0.05*(max(to.plot$x)), max.mutation.time)) +
      ggtitle(i)
    
    print(p)
    
    p <- ggplot(to.plot[to.plot$Timing=="Pre-CNV",], aes(x=x, fill=Segment)) + geom_histogram( binwidth = binwidth)+ 
      scale_fill_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot[to.plot$Timing=="Pre-CNV",])))+
      geom_density(data=to.plot[to.plot$Timing=="Post-CNV",], linetype=2, aes(x=x), inherit.aes = F) +
      scale_y_continuous(name="# Genomic segments") + 
        scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$x)), max(to.plot$x)*1.1)) +
      ggtitle(i)
    
    print(p)
    
    to.plot. <- to.plot[to.plot$Timing=="Pre-CNV",]
    to.plot.$Segment <- factor(to.plot.$Segment, levels=unique(to.plot.$Segment))
    
    p <- ggplot(to.plot., aes(x=x, xmin=xlower, xmax=xupper, y=as.numeric(Segment)/nrow(to.plot.), col=Segment)) +
      scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot.))) + 
      geom_point() + geom_errorbarh(height=0)+ 
      geom_vline(data=data.frame(x=mutation.time.mrca[i]), aes(xintercept=x/3.3/10^3, y=1), inherit.aes = F, col=manual.colors["Late"],
                 linetype=2) + 
      geom_ribbon(data=data.frame(xmin=rep(mutation.time.mrca.lower[i],2),
                                  xmax=rep(mutation.time.mrca.upper[i],2),
                                  y=c(0, nrow(to.plot[to.plot$Timing=="Pre-CNV",]))), 
                  aes(xmin=xmin/3.3/10^3, xmax=xmax/3.3/10^3,y=y),
                  inherit.aes = F, fill=manual.colors["Late"], col=NA, alpha=0.5) +
           scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$x)), max(to.plot$xupper)*1.1)) +
      scale_y_continuous(limits=c(0, nrow(to.plot)))
    
    print(p)   
    
    dev.off()
  }
  
  if(i == "NBE99"){
    pdf(paste0(panel.directory, "Figure_S2c.pdf"))
    
    p <- ggplot(to.plot[to.plot$Timing=="Post-CNV",], aes(x=x)) + geom_histogram( binwidth = binwidth, fill=manual.colors["Late"])+ 
      scale_y_continuous(name="# Genomic segments") + 
         scale_x_continuous(name="Mutations per Mb", limits = c(-0.05*(max(to.plot$x)), max.mutation.time)) +
      ggtitle(i)
    
    print(p)
    
    p <- ggplot(to.plot[to.plot$Timing=="Pre-CNV",], aes(x=x, fill=Segment)) + geom_histogram( binwidth = binwidth)+ 
      scale_fill_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot[to.plot$Timing=="Pre-CNV",])))+
      geom_density(data=to.plot[to.plot$Timing=="Post-CNV",], linetype=2, aes(x=x), inherit.aes = F) +
      scale_y_continuous(name="# Genomic segments") + 
      scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$x)), max(to.plot$x)*1.1)) +
      ggtitle(i)
    
    print(p)
    
    to.plot. <- to.plot[to.plot$Timing=="Pre-CNV",]
    to.plot.$Segment <- factor(to.plot.$Segment, levels=unique(to.plot.$Segment))
    
    p <- ggplot(to.plot., aes(x=x, xmin=xlower, xmax=xupper, y=as.numeric(Segment)/nrow(to.plot.), col=Segment)) +
      scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot.))) + 
      geom_point() + geom_errorbarh(height=0)+ 
      geom_vline(data=data.frame(x=mutation.time.mrca[i]), aes(xintercept=x/3.3/10^3, y=1), inherit.aes = F, col=manual.colors["Late"],
                 linetype=2) + 
      geom_ribbon(data=data.frame(xmin=rep(mutation.time.mrca.lower[i],2),
                                  xmax=rep(mutation.time.mrca.upper[i],2),
                                  y=c(0, nrow(to.plot[to.plot$Timing=="Pre-CNV",]))), 
                  aes(xmin=xmin/3.3/10^3, xmax=xmax/3.3/10^3,y=y),
                  inherit.aes = F, fill=manual.colors["Late"], col=NA, alpha=0.5) +
      scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$x)), max(to.plot$xupper)*1.1)) +
      scale_y_continuous(limits=c(0, nrow(to.plot)))
    
    print(p)    

    dev.off()
  }
  
}



##########################################################################################################################################
#### Figure S2a: # Mutation densities per segment for tumor NBE15

example.tumor <- "NBE15"

to.plot <- data.frame(x= mutation.time.most.likely[[example.tumor]]/3.3/10^3,
                      xupper= mutation.time.upper[[example.tumor]]/3.3/10^3,
                      xlower=mutation.time.lower[[example.tumor]]/3.3/10^3,
                      Segment = names(mutation.time.most.likely[[example.tumor]]),
                      Timing = sapply(names(mutation.time.most.likely[[example.tumor]]), function(x){
                        x <- strsplit(x, split="_")[[1]]
                        if(is.na(x[4])){
                          return("Post-CNV")
                        }
                        if(x[4]=="I"){
                          return("Post-CNV")
                        }else{
                          return("Pre-CNV")
                        }
                      }))
to.plot$Segment <- factor(to.plot$Segment, levels=unique(to.plot$Segment))

to.plot. <- to.plot[to.plot$Timing=="Post-CNV",]
to.plot.$Segment <- factor(to.plot.$Segment, levels=unique(to.plot.$Segment))
to.plot.$CopyNumber <- sapply(to.plot.$Segment, function(x){
  strsplit(as.character(x), split="_")[[1]][3]
})


to.plot.$CopyNumber <- factor(to.plot.$CopyNumber, levels=c("haploid", "disomic", "trisomic", "tetrasomic"))


addWorksheet(wb.s, "a")
writeData(wb.s, "a", to.plot.)


pdf(paste0(panel.directory, "Figure_S2a.pdf"), width=3, height=3)

p <- ggplot(to.plot., aes(x=CopyNumber, y=x, ymin=xlower, ymax=xupper)) + geom_boxplot() +
  geom_beeswarm()+
  #geom_pointrange(position = position_jitter(w = 0.1, h = 0)) + 
  scale_y_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$x)), max(to.plot$xupper)*1.1)) +
  geom_signif(comparisons=list(c("disomic", "trisomic"), c("disomic", "tetrasomic"), c("trisomic", "tetrasomic"))) +
  ggtitle(sample.information.80x[i,]$Evolution_paper_Id) + scale_y_continuous(name="")

print(p)
#ggarrange(plotlist=p.mrca.densities.all.tumors, nrow=6, ncol=4)

dev.off()


##############################################################################################################################################
## Figure S2b_d_f_h: density distributions for MRCA and ECA - detection data set
source(paste0(custom.script.directory, "Oncoprint.R"))

## if ECA was not uniquely identified, take the earliest time point
mutation.time.eca[names(earliest.mutation.time)] <- earliest.mutation.time
mutation.time.eca.lower[names(earliest.mutation.time)] <- earliest.mutation.time.lower
mutation.time.eca.upper[names(earliest.mutation.time)] <- earliest.mutation.time.upper
####

colnames(sample.information.80x)[which(colnames(sample.information.80x)=="ECA")] <- "ECA.exists"

max.mutation.time.primary <- max(mutation.time.mrca[sample.information.80x$Patient_ID[sample.information.80x$Location=="Primary" &
                                                                                               sample.information.80x$Treatment==FALSE]], na.rm=T)

sample.information.80x$Telomere.maintenance.mechanism <- factor(sample.information.80x$Telomere.maintenance.mechanism,
                                                                levels=c("MNA", "TERT", "ALT", "Multiple", "None"))

##### Primary tumors, with and without treatment, metastases
## Plot different ploidies with ECA, merge di- and tetraploids

p1 <- list()
p1.subset.list <- list()

for(ploidy in c("triploid", "di-tetraploid")){
  
  if(ploidy=="triploid"){
    ploidy <- 3
  }else{
    ploidy <- c(2,4)
  }
  
  subset=sample.information.80x[ sample.information.80x$Location %in% c("Primary", "Metastasis") &
                                   sample.information.80x$ECA.exists==T &
                                   sample.information.80x$Ploidy %in% ploidy,,drop=F]
  
  if(nrow(subset)==0){next}
  
  
  subset$Telomere.maintenance.mechanism <- factor(subset$Telomere.maintenance.mechanism,
                                                  levels=c("MNA", "TERT", "ALT", "Multiple", "None"))
  
  subset$MRCAtime <- mutation.time.mrca[subset$Patient_ID]
  subset <- subset[order(subset$MRCAtime),]
  
  
  to.plot <- cbind(subset, data.frame(MRCA=mutation.time.mrca[subset$Patient_ID]/3.3/10^3,
                                      ECA=mutation.time.eca[subset$Patient_ID]/3.3/10^3,
                                      MRCA.upper=mutation.time.mrca.upper[subset$Patient_ID]/3.3/10^3,
                                      MRCA.lower=mutation.time.mrca.lower[subset$Patient_ID]/3.3/10^3,
                                      ECA.upper=mutation.time.eca.upper[subset$Patient_ID]/3.3/10^3,
                                      ECA.lower=mutation.time.eca.lower[subset$Patient_ID]/3.3/10^3))
  
  
  
  if(ploidy%in%c(2,4)){
    panel="d_population_summary"
  }else if(ploidy==3){
    panel="b"
  }
  
  addWorksheet(wb.s, panel)
  writeData(wb.s, panel, to.plot[,c("ECA", "ECA.lower", "ECA.upper", "MRCA", "MRCA.lower", "MRCA.upper")])
  
  
  
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
    scale_x_continuous(name = "Mutations/Mb",
                       limits=c(0, max.mutation.time.primary/3.3/10^3))+
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(name = "Fraction of tumors") +
    ggtitle(paste("# Cases = ", nrow(subset))) + 
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
  
  p1.subset.list[[length(p1.subset.list) + 1]] <- ggplot(data=to.plot.,
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
    
    p1.subset.list[[length(p1.subset.list)]] <- p1.subset.list[[length(p1.subset.list)]] +
      geom_rect(data=to.plot., aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill=chr.change), inherit.aes=T)
  }
  
  
  
}


## Plot the cases without ECA separately per ploidy; merge di- and tetraploids
p2 <- list()
p2.subset.list <- list()

for(ploidy in c("triploid", "di-tetraploid")){
  
  if(ploidy=="triploid"){
    ploidy <- 3
  }else{
    ploidy <- c(2,4)
  }
  
  subset=sample.information.80x[ sample.information.80x$Location %in% c("Primary", "Metastasis") &
                                   sample.information.80x$ECA.exists==F &
                                   sample.information.80x$Ploidy %in% ploidy,,drop=F]
  
  if(nrow(subset)==0){next}
  
  
  subset$Telomere.maintenance.mechanism <- factor(subset$Telomere.maintenance.mechanism,
                                                  levels=c("MNA", "TERT", "ALT", "Multiple", "None"))
  
  subset$MRCAtime <- mutation.time.mrca[subset$Patient_ID]
  subset <- subset[order(subset$MRCAtime),]
  
  
  to.plot <- cbind(subset, data.frame(MRCA=mutation.time.mrca[subset$Patient_ID]/3.3/10^3,
                                      ECA=mutation.time.eca[subset$Patient_ID]/3.3/10^3,
                                      MRCA.upper=mutation.time.mrca.upper[subset$Patient_ID]/3.3/10^3,
                                      MRCA.lower=mutation.time.mrca.lower[subset$Patient_ID]/3.3/10^3,
                                      ECA.upper=mutation.time.eca.upper[subset$Patient_ID]/3.3/10^3,
                                      ECA.lower=mutation.time.eca.lower[subset$Patient_ID]/3.3/10^3))
  
  
  
  if(ploidy%in%c(2,4)){
    panel="e_population_summary"
  }else if(ploidy==3){
    panel="c_population_summary"
  }
  
  addWorksheet(wb.s, panel)
  writeData(wb.s, panel, to.plot[,c("ECA", "ECA.lower", "ECA.upper", "MRCA", "MRCA.lower", "MRCA.upper")])
  
  
  
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
  
  p2.subset.list[[length(p2.subset.list) + 1]] <- ggplot(data=to.plot.,
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
    
    p2.subset.list[[length(p2.subset.list)]] <- p2.subset.list[[length(p2.subset.list)]] +
      geom_rect(data=to.plot., aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill=chr.change), inherit.aes=T)
  }
  
  
  
}


pdf(paste0(panel.directory, "Figure_S2_b_c_d_e_population_summary.pdf"), width=9, height=9, useDingbats = F)

figure <- ggarrange(plotlist=c(p1, p2), nrow=3, ncol=3)
annotate_figure(figure, top="Primary tumor / Metastasis")

figure <- ggarrange(plotlist=c(p1.subset.list, p2.subset.list), nrow=3, ncol=3)
annotate_figure(figure, top="Primary tumor / Metastasis")

dev.off()


##########################################################################################################################################
#### Figure S2f: # of events for triploidization

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


pdf(paste0(panel.directory,"Figure_S2f.pdf"), width=2, height=2.2, useDingbats = F)

to.plot <- rbind(data.frame(Chromosomes=numbers.of.3n.chromosomes.at.single.event/numbers.of.3n.chromosomes,
                            Tumor=triploid.tumors.80x,
                            Telomere.type=telomere.classification.80x[triploid.tumors.80x],
                            Time = sample.information.80x[triploid.tumors.80x,]$Location,
                            Ploidy=3),
                 data.frame(Chromosomes=numbers.of.4n.chromosomes.at.single.event/numbers.of.4n.chromosomes,
                            Tumor=tetraploid.tumors.80x,
                            Telomere.type=telomere.classification.80x[tetraploid.tumors.80x],
                            Time = sample.information.80x[tetraploid.tumors.80x,]$Location,
                            Ploidy=4))


panel="f"
addWorksheet(wb.s, panel)
writeData(wb.s, panel, to.plot)


ggplot(to.plot, aes(x=Ploidy, y=Chromosomes,group=Ploidy)) + geom_boxplot(width=0.5) + geom_beeswarm() + scale_y_continuous(limits=c(0,1.01))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


dev.off()



##############################################################################################################################################
##### Figure S2g: Compare ECA and MRCA between primary and relapse tumor

mutation.time.eca[names(earliest.mutation.time)] <- earliest.mutation.time
mutation.time.eca.lower[names(earliest.mutation.time)] <- earliest.mutation.time.lower
mutation.time.eca.upper[names(earliest.mutation.time)] <- earliest.mutation.time.upper

colnames(sample.information.80x)[which(colnames(sample.information.80x)=="ECA")] <- "ECA.exists"
to.plot <- cbind(sample.information.80x, data.frame(MRCA=mutation.time.mrca[rownames(sample.information.80x)]/3.3/10^3,
                                                    ECA=mutation.time.eca[rownames(sample.information.80x)]/3.3/10^3,
                                                    MRCA.upper=mutation.time.mrca.upper[rownames(sample.information.80x)]/3.3/10^3,
                                                    MRCA.lower=mutation.time.mrca.lower[rownames(sample.information.80x)]/3.3/10^3,
                                                    ECA.upper=mutation.time.eca.upper[rownames(sample.information.80x)]/3.3/10^3,
                                                    ECA.lower=mutation.time.eca.lower[rownames(sample.information.80x)]/3.3/10^3))

to.plot$Time <- sample.information.80x$Location
to.plot$Time[to.plot$Time=="Relapse 3"] <- "Relapse"
to.plot$Time[to.plot$Time %in% c("Relapse tumor", "Relapse metastasis")] <- "Relapse"

to.plot$Subtype <- sample.information.80x$Subtype
to.plot$Time <- factor(to.plot$Time, levels=c("Primary", "Metastasis", "Relapse"))


addWorksheet(wb.s, "g_ECA")
writeData(wb.s, "g_ECA", to.plot[,c("ECA", "ECA.lower", "ECA.upper", "Ploidy", "Time")])
addWorksheet(wb.s, "g_MRCA")
writeData(wb.s, "g_MRCA", to.plot[,c("MRCA", "MRCA.lower", "MRCA.upper", "Ploidy", "Time")])


pdf(paste0(panel.directory, "Figure_S2g.pdf"), width=4, height=4, useDingbats = F)

ggplot(to.plot[to.plot$ECA!=to.plot$MRCA,], aes(x=Time, y=ECA, col=as.character(Ploidy)))+geom_quasirandom() +
  scale_color_manual(values=c("2" = "black", "3" = "orange", "4" = "firebrick")) +
  scale_x_discrete( breaks = c("Primary FALSE", "Primary TRUE", "Metastasis FALSE"),
                    labels=c("Primary untreated", "Primary treated", "Metastasis"))+ 
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1,
         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_y_continuous(name = "SSNVs/Mb at ECA")

wilcox.test(to.plot[to.plot$Time=="Primary" & to.plot$ECA!= to.plot$MRCA,]$ECA, to.plot[to.plot$Time=="Metastasis" & to.plot$ECA!= to.plot$MRCA,]$ECA)
wilcox.test(to.plot[to.plot$Time=="Primary" & to.plot$ECA!= to.plot$MRCA,]$ECA, to.plot[to.plot$Time=="Relapse" & to.plot$ECA!= to.plot$MRCA,]$ECA)
wilcox.test(to.plot[to.plot$Time=="Relapse" & to.plot$ECA!= to.plot$MRCA,]$ECA, to.plot[to.plot$Time=="Metastasis" & to.plot$ECA!= to.plot$MRCA,]$ECA)

ggplot(to.plot, aes(x=Time, y=MRCA, col=as.character(Ploidy)))+geom_quasirandom() +
  scale_color_manual(values=c("2" = "black", "3" = "orange", "4" = "firebrick")) +
  scale_x_discrete( breaks = c("Primary FALSE", "Primary TRUE", "Metastasis FALSE"),
                    labels=c("Primary untreated", "Primary treated", "Metastasis"))+ 
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1,
         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_y_continuous(name = "SSNVs/Mb at MRCA")

wilcox.test(to.plot[to.plot$Time=="Primary",]$MRCA, to.plot[to.plot$Time=="Metastasis",]$MRCA)
wilcox.test(to.plot[to.plot$Time=="Primary",]$MRCA, to.plot[to.plot$Time=="Relapse",]$MRCA)
wilcox.test(to.plot[to.plot$Time=="Metastasis",]$MRCA, to.plot[to.plot$Time=="Relapse",]$MRCA)

dev.off()



##########################################################################################################################################
saveWorkbook(wb, file = paste0(panel.directory, "Source_data_Fig.2.xlsx"), overwrite=T)

saveWorkbook(wb.s, file = paste0(panel.directory, "Source_data_Fig.S2.xlsx"), overwrite=T)



##############################################################################################################################################
## Figure S3: density distributions for MRCA and ECA - validation data set

## source data:
wb.s <- createWorkbook()

source(paste0(custom.script.directory, "Oncoprint_additional_samples.R"))

## if ECA was not uniquely identified, take the earliest time point
mutation.time.eca[names(earliest.mutation.time)] <- earliest.mutation.time
mutation.time.eca.lower[names(earliest.mutation.time)] <- earliest.mutation.time.lower
mutation.time.eca.upper[names(earliest.mutation.time)] <- earliest.mutation.time.upper
####

colnames(sample.information.30x)[which(colnames(sample.information.30x)=="ECA")] <- "ECA.exists"
sample.information.30x$ECA.exists <- as.logical(sample.information.30x$ECA.exists)
sample.information.30x$Ploidy <- sample.information.30x$Rounded.ploidy

max.mutation.time.primary <- max(mutation.time.mrca[rownames(sample.information.30x)[sample.information.30x$Location=="Primary"]], na.rm=T)

sample.information.30x$Telomere.maintenance.mechanism <- factor(sample.information.30x$Telomere.maintenance.mechanism,
                                                                levels=c("MNA", "TERT", "ALT", "Multiple", "None"))

##### Primary tumors, with and without treatment, metastases
## Plot different ploidies with ECA, merge di- and tetraploids

p1 <- list()
p1.subset.list <- list()

for(ploidy in c("triploid", "di-tetraploid")){
  
  if(ploidy=="triploid"){
    ploidy <- 3
  }else{
    ploidy <- c(2,4)
  }
  
  subset=sample.information.30x[ sample.information.30x$Location %in% c("Primary", "Metastasis") &
                                   sample.information.30x$ECA.exists==T &
                                   sample.information.30x$Ploidy %in% ploidy,,drop=F]
  
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
  
  
  
  if(ploidy%in%c(2,4)){
    panel="Di_tetraploid_ECA"
  }else if(ploidy==3){
    panel="Triploid_ECA"
  }
  
  addWorksheet(wb.s, panel)
  writeData(wb.s, panel, to.plot[,c("ECA", "ECA.lower", "ECA.upper", "MRCA", "MRCA.lower", "MRCA.upper")])
  
  
  
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
    scale_x_continuous(name = "Mutations/Mb",
                       limits=c(0, max.mutation.time.primary/3.3/10^3))+
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(name = "Fraction of tumors") +
    ggtitle(paste("# Cases = ", nrow(subset))) + 
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
  
  
  ## TMM, Clinical subtype & CNVs
  to.plot. <- subset[,c("MRCAtime", "Telomere.maintenance.mechanism", "ManualScore", "Ploidy")]
  to.plot.$xmin <-c(0, (to.plot.$MRCAtime/3.3/10^3)[-length(to.plot.$MRCAtime)])
  to.plot.$xmax <- to.plot.$MRCAtime/3.3/10^3
  to.plot.$ymin <- 0
  to.plot.$ymax <- 0.23
  
  to.plot.[,c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2")] <- t(matrix.early.late.30x[c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2"),
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
  
  p1.subset.list[[length(p1.subset.list) + 1]] <- ggplot(data=to.plot.,
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
    
    p1.subset.list[[length(p1.subset.list)]] <- p1.subset.list[[length(p1.subset.list)]] +
      geom_rect(data=to.plot., aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill=chr.change), inherit.aes=T)
  }
  
  
  
}


## Plot the cases without ECA separately per ploidy; merge di- and tetraploids
p2 <- list()
p2.subset.list <- list()

for(ploidy in c("triploid", "di-tetraploid")){
  
  if(ploidy=="triploid"){
    ploidy <- 3
  }else{
    ploidy <- c(2,4)
  }
  
  subset=sample.information.30x[ sample.information.30x$Location %in% c("Primary", "Metastasis") &
                                   sample.information.30x$ECA.exists==F &
                                   sample.information.30x$Ploidy %in% ploidy,,drop=F]
  
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
  
  
  
  if(ploidy%in%c(2,4)){
    panel="Di_tetraploid_no_ECA"
  }else if(ploidy==3){
    panel="Triploid_no_ECA"
  }
  
  addWorksheet(wb.s, panel)
  writeData(wb.s, panel, to.plot[,c("ECA", "ECA.lower", "ECA.upper", "MRCA", "MRCA.lower", "MRCA.upper")])
  
  
  
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
  
  to.plot.[,c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2")] <- t(matrix.early.late.30x[c("17q", "1p", "1q", "7q", "2p", "11q", "17", "1", "7", "2"),
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
  
  p2.subset.list[[length(p2.subset.list) + 1]] <- ggplot(data=to.plot.,
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
    
    p2.subset.list[[length(p2.subset.list)]] <- p2.subset.list[[length(p2.subset.list)]] +
      geom_rect(data=to.plot., aes_string(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill=chr.change), inherit.aes=T)
  }
  
  
  
}


pdf(paste0(panel.directory, "Figure_S3.pdf"), width=9, height=9, useDingbats = F)

figure <- ggarrange(plotlist=c(p1, p2), nrow=3, ncol=3)
annotate_figure(figure, top="Primary tumor / Metastasis")

figure <- ggarrange(plotlist=c(p1.subset.list, p2.subset.list), nrow=3, ncol=3)
annotate_figure(figure, top="Primary tumor / Metastasis")

dev.off()


##############################################################################################################################################
saveWorkbook(wb.s, file = paste0(panel.directory,"Source_data_Fig.3.xlsx"), overwrite=T)


