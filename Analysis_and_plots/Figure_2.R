## Reproduce Fig. 2
##############################################################################################################################################
## load settings and libraries

source("./Settings.R")
## either re-run
#source(paste0(custom.script.directory, "Mutation_density_quantification.R"))
## or load
load(paste0(rdata.directory, "Purity_ploidy.RData"))
load(paste0(rdata.directory, "Clonal_mutations_different_CNs.RData"))
load(paste0(rdata.directory, "MRCA_timing.RData"))
load(paste0(rdata.directory, "Vafs_all_tumors.RData"))

## get mutation rate per day from later analysis (Figure 5; code this differently/in order)

## source data:
wb <- createWorkbook()
wb.s2 <- createWorkbook()
wb.s3 <- createWorkbook()

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

##########################################################################################################################################
#### Figure S2a: # Mutation densities per segment for tumor NBE15

example.tumor <- "NBE15"

to.plot <- mutation.time[[example.tumor]]

to.plot$Timing <- sapply(as.character(to.plot$Segment), function(x){
                        x <- strsplit(x, split="_")[[1]]
                        if(is.na(x[4])){
                          return("Post-CNV")
                        }
                        if(x[4]=="I"){
                          return("Post-CNV")
                        }else{
                          return("Pre-CNV")
                        }
                      })
to.plot$Segment <- factor(to.plot$Segment, levels=unique(to.plot$Segment))

to.plot. <- to.plot[to.plot$Timing=="Post-CNV",]
to.plot.$Segment <- factor(to.plot.$Segment, levels=unique(to.plot.$Segment))
to.plot.$CopyNumber <- sapply(to.plot.$Segment, function(x){
  strsplit(as.character(x), split="_")[[1]][3]
})


to.plot.$CopyNumber <- factor(to.plot.$CopyNumber, levels=c("haploid", "disomic", "trisomic", "tetrasomic"))

to.plot.$Mean <- to.plot.$Mean/3.3/10^3
to.plot.$Max <- to.plot.$Max/3.3/10^3
to.plot.$Min <- to.plot.$Min/3.3/10^3

addWorksheet(wb.s2, "a")
writeData(wb.s2, "a", to.plot.)


pdf(paste0(panel.directory, "Figure_S2a.pdf"), width=3, height=3)

p <- ggplot(to.plot., aes(x=CopyNumber, y=Mean, ymin=Min, ymax=Max)) + geom_boxplot() +
  geom_beeswarm()+
  scale_y_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$Mean)), max(to.plot$Max)*1.1)) +
  geom_signif(comparisons=list(c("disomic", "trisomic"), c("disomic", "tetrasomic"), c("trisomic", "tetrasomic"))) +
  ggtitle(sample.information.discovery[example.tumor,]$Evolution_paper_Id) + scale_y_continuous(name="")

print(p)

## do tests:

chars <- capture.output(wilcox.test(to.plot.$Mean[to.plot.$CopyNumber=="disomic"],
                                    to.plot.$Mean[to.plot.$CopyNumber=="trisomic"], conf.int = T))

writeData(wb.s2, sheet = "a", chars, startCol = ncol(to.plot.)+5)

chars <- capture.output(wilcox.test(to.plot.$Mean[to.plot.$CopyNumber=="disomic"], 
                                    to.plot.$Mean[to.plot.$CopyNumber=="tetrasomic"], conf.int = T))

writeData(wb.s2, sheet = "a", chars, startCol = ncol(to.plot.)+7)

chars <- capture.output(wilcox.test(to.plot.$Mean[to.plot.$CopyNumber=="trisomic"], 
                                    to.plot.$Mean[to.plot.$CopyNumber=="tetrasomic"], conf.int = T))

writeData(wb.s2, sheet = "a", chars, startCol = ncol(to.plot.)+9)

dev.off()

##############################################################################################################################################
## Figure 2d-f and Fig. S3 a, c, d: density plots per genomic segment for NBE15, NBE14, NBE44, NBE99

## manually set the range 
max.mutation.time <- 0.7

for(i in c("NBE15", "NBE14", "NBE44", "NBE99")){

  binwidth = (max(mutation.time[[i]]$Mean/3.3/10^3) - min(mutation.time[[i]]$Mean/3.3/10^3))/20
  
  to.plot <- mutation.time[[i]]
  to.plot$Mean <- to.plot$Mean/3.3/10^3
  to.plot$Max <- to.plot$Max/3.3/10^3
  to.plot$Min <- to.plot$Min/3.3/10^3
  to.plot$Timing <- sapply(as.character(to.plot$Segment), function(x){
                          x <- strsplit(x, split="_")[[1]]
                          if(is.na(x[4])){
                            return("Post-CNV")
                          }
                          if(x[4]=="I"){
                            return("Post-CNV")
                          }else{
                            return("Pre-CNV")
                          }
                        })
  
  to.plot$Segment <- factor(to.plot$Segment, levels=unique(to.plot$Segment))
  ## for the three examples shown in the manuscript, export the source data
  if(i %in% c("NBE15", "NBE14", "NBE44", "NBE99")){
    
    if(i=="NBE15"){
      panel="d_e_f"
      addWorksheet(wb, panel)
      writeData(wb, panel, to.plot)
      addWorksheet(wb, "f_ECA_MRCA")
      source.data.eca.mrca <- data.frame(ECA = mutation.time.eca[i,]$Mean/3.3/10^3,
                                         ECA.lower = mutation.time.eca[i,]$Min/3.3/10^3,
                                         ECA.higher = mutation.time.eca[i,]$Max/3.3/10^3,
                                         MRCA = mutation.time.mrca[i,]$Mean/3.3/10^3,
                                         MRCA.lower = mutation.time.mrca[i,]$Min/3.3/10^3,
                                         MRCA.higher = mutation.time.mrca[i,]$Max/3.3/10^3)
      writeData(wb, "f_ECA_MRCA", source.data.eca.mrca)
      
      p.values <- mrca.eca[[i]]$p.values.mrca
      p.values$Chromosome <- sapply(rownames(p.values), function(x){
        strsplit(x, split="_")[[1]][2]
      })
      p.values$CN  <- sapply(rownames(p.values), function(x){
        x <- strsplit(x, split="_")[[1]][3]
        if(x=="monosomic"){
          1
        }else if(x=="disomic"){
          2
        }else if(x=="trisomic"){
          3
        }else if(x=="tetrasomic"){
          4
        }
      })
      p.values$A  <- sapply(rownames(p.values), function(x){
        strsplit(x, split="_")[[1]][4]
      })
      ## restrict to gains
      p.values <- p.values[p.values$A!="I",c("Chromosome", "CN", "p", "adj.p")]
      p.values$sig <- ifelse(p.values$adj.p < 0.01, "**", "n.s.")
      rownames(p.values) <- NULL
      
      writeData(wb, sheet = "f_ECA_MRCA", p.values, startCol = ncol(source.data.eca.mrca)+5, rowNames = F)
      
      
    }else if(i=="NBE14"){
      panel="c"
      addWorksheet(wb.s3, panel)
      writeData(wb.s3, panel, to.plot)
      
      addWorksheet(wb.s3, "c_ECA_MRCA")
      source.data.eca.mrca <- data.frame(ECA = mutation.time.eca[i,]$Mean/3.3/10^3,
                                         ECA.lower = mutation.time.eca[i,]$Min/3.3/10^3,
                                         ECA.higher = mutation.time.eca[i,]$Max/3.3/10^3,
                                         MRCA = mutation.time.mrca[i,]$Mean/3.3/10^3,
                                         MRCA.lower = mutation.time.mrca[i,]$Min/3.3/10^3,
                                         MRCA.higher = mutation.time.mrca[i,]$Max/3.3/10^3)
      writeData(wb.s3, "c_ECA_MRCA", source.data.eca.mrca)
      
      p.values <- mrca.eca[[i]]$p.values.mrca
      p.values$Chromosome <- sapply(rownames(p.values), function(x){
        strsplit(x, split="_")[[1]][2]
      })
      p.values$CN  <- sapply(rownames(p.values), function(x){
        x <- strsplit(x, split="_")[[1]][3]
        if(x=="monosomic"){
          1
        }else if(x=="disomic"){
          2
        }else if(x=="trisomic"){
          3
        }else if(x=="tetrasomic"){
          4
        }
      })
      p.values$A  <- sapply(rownames(p.values), function(x){
        strsplit(x, split="_")[[1]][4]
      })
      ## restrict to gains
      p.values <- p.values[p.values$A!="I",c("Chromosome", "CN", "p", "adj.p")]
      p.values$sig <- ifelse(p.values$adj.p < 0.01, "**", "n.s.")
      rownames(p.values) <- NULL
      
      writeData(wb.s3, sheet = "c_ECA_MRCA", p.values, startCol = ncol(source.data.eca.mrca)+5, rowNames = F)
    }else if(i=="NBE99"){
      panel="a"
      addWorksheet(wb.s3, panel)
      writeData(wb.s3, panel, to.plot)
      
      addWorksheet(wb.s3, "a_ECA_MRCA")
      source.data.eca.mrca <- data.frame(ECA = mutation.time.eca[i,]$Mean/3.3/10^3,
                                         ECA.lower = mutation.time.eca[i,]$Min/3.3/10^3,
                                         ECA.higher = mutation.time.eca[i,]$Max/3.3/10^3,
                                         MRCA = mutation.time.mrca[i,]$Mean/3.3/10^3,
                                         MRCA.lower = mutation.time.mrca[i,]$Min/3.3/10^3,
                                         MRCA.higher = mutation.time.mrca[i,]$Max/3.3/10^3)
      writeData(wb.s3, "a_ECA_MRCA", source.data.eca.mrca)
      
      p.values <- mrca.eca[[i]]$p.values.mrca
      p.values$Chromosome <- sapply(rownames(p.values), function(x){
        strsplit(x, split="_")[[1]][2]
      })
      p.values$CN  <- sapply(rownames(p.values), function(x){
        x <- strsplit(x, split="_")[[1]][3]
        if(x=="monosomic"){
          1
        }else if(x=="disomic"){
          2
        }else if(x=="trisomic"){
          3
        }else if(x=="tetrasomic"){
          4
        }
      })
      p.values$A  <- sapply(rownames(p.values), function(x){
        strsplit(x, split="_")[[1]][4]
      })
      ## restrict to gains
      p.values <- p.values[p.values$A!="I",c("Chromosome", "CN", "p", "adj.p")]
      p.values$sig <- ifelse(p.values$adj.p < 0.01, "**", "n.s.")
      rownames(p.values) <- NULL
      
      writeData(wb.s3, sheet = "a_ECA_MRCA", p.values, startCol = ncol(source.data.eca.mrca)+5, rowNames = F)
    }else{
      panel="d"
      addWorksheet(wb.s3, panel)
      writeData(wb.s3, panel, to.plot)
      
      addWorksheet(wb.s3, "d_ECA_MRCA")
      source.data.eca.mrca <- data.frame(ECA = mutation.time.eca[i,]$Mean/3.3/10^3,
                                         ECA.lower = mutation.time.eca[i,]$Min/3.3/10^3,
                                         ECA.higher = mutation.time.eca[i,]$Max/3.3/10^3,
                                         MRCA = mutation.time.mrca[i,]$Mean/3.3/10^3,
                                         MRCA.lower = mutation.time.mrca[i,]$Min/3.3/10^3,
                                         MRCA.higher = mutation.time.mrca[i,]$Max/3.3/10^3)
      writeData(wb.s3, "d_ECA_MRCA", source.data.eca.mrca)
      
      p.values <- mrca.eca[[i]]$p.values.mrca
      p.values$Chromosome <- sapply(rownames(p.values), function(x){
        strsplit(x, split="_")[[1]][2]
      })
      p.values$CN  <- sapply(rownames(p.values), function(x){
        x <- strsplit(x, split="_")[[1]][3]
        if(x=="monosomic"){
          1
        }else if(x=="disomic"){
          2
        }else if(x=="trisomic"){
          3
        }else if(x=="tetrasomic"){
          4
        }
      })
      p.values$A  <- sapply(rownames(p.values), function(x){
        strsplit(x, split="_")[[1]][4]
      })
      ## restrict to gains
      p.values <- p.values[p.values$A!="I",c("Chromosome", "CN", "p", "adj.p")]
      p.values$sig <- ifelse(p.values$adj.p < 0.01, "**", "n.s.")
      rownames(p.values) <- NULL
      
      writeData(wb.s3, sheet = "d_ECA_MRCA", p.values, startCol = ncol(source.data.eca.mrca)+5, rowNames = F)
    }
    

    
  }
 
  if(i == "NBE15"){
    pdf(paste0(panel.directory, "Figure_2d_e_f.pdf"))
    
    p <- ggplot(to.plot[to.plot$Timing=="Post-CNV",], aes(x=Mean)) + geom_histogram( binwidth = binwidth, fill=manual.colors["Late"])+ 
      scale_y_continuous(name="# Genomic segments") + 
       scale_x_continuous(name="Mutations per Mb", limits = c(-0.05*(max(to.plot$Mean)), max.mutation.time)) +
      ggtitle(i)
    
    print(p)
    
    p <- ggplot(to.plot[to.plot$Timing=="Pre-CNV",], aes(x=Mean, fill=Segment)) + geom_histogram( binwidth = binwidth)+ 
      scale_fill_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot[to.plot$Timing=="Pre-CNV",])))+
      geom_density(data=to.plot[to.plot$Timing=="Post-CNV",], linetype=2, aes(x=Mean), inherit.aes = F) +
      scale_y_continuous(name="# Genomic segments") + 
        scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$Mean)), max(to.plot$Mean)*1.1)) +
      ggtitle(i)
    
    print(p)
    
    to.plot. <- to.plot[to.plot$Timing=="Pre-CNV",]
    to.plot.$Segment <- factor(to.plot.$Segment, levels=unique(to.plot.$Segment))
    
    p <- ggplot(to.plot., aes(x=Mean, xmin=Min, xmax=Max, y=as.numeric(Segment)/nrow(to.plot.), col=Segment)) +
      scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot.))) + 
      geom_point() + geom_errorbarh(height=0)+ 
      geom_vline(data=data.frame(x=mutation.time.mrca[i,]$Mean), aes(xintercept=x/3.3/10^3, y=1), inherit.aes = F, col=manual.colors["Late"],
                 linetype=2) + 
      geom_ribbon(data=data.frame(xmin=rep(mutation.time.mrca[i,]$Min,2),
                                  xmax=rep(mutation.time.mrca[i,]$Max,2),
                                  y=c(0, nrow(to.plot[to.plot$Timing=="Pre-CNV",]))), 
                  aes(xmin=xmin/3.3/10^3, xmax=xmax/3.3/10^3,y=y),
                  inherit.aes = F, fill=manual.colors["Late"], col=NA, alpha=0.5) +
      geom_vline(data=data.frame(x=mutation.time.eca[i,]$Mean), aes(xintercept=x/3.3/10^3, y=1), inherit.aes = F, col=manual.colors["Early"],
                 linetype=2) + 
      geom_ribbon(data=data.frame(xmin=rep(mutation.time.eca[i,]$Min, 2),
                                  xmax=rep(mutation.time.eca[i,]$Max,2),
                                  y=c(0, nrow(to.plot[to.plot$Timing=="Pre-CNV",]))), 
                  aes( xmin=xmin/3.3/10^3, xmax=xmax/3.3/10^3, y=y),
                  height=0, inherit.aes = F, fill=manual.colors["Early"], col=NA, alpha=0.5) +
      scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$Mean)), max(to.plot$Max)*1.1)) +
      scale_y_continuous(limits=c(0, nrow(to.plot)))
    
    print(p)

    dev.off()
    
  }
  
  if(i == "NBE14"){
    pdf(paste0(panel.directory, "Figure_S3c.pdf"))
    
    p <- ggplot(to.plot[to.plot$Timing=="Post-CNV",], aes(x=Mean)) + geom_histogram( binwidth = binwidth, fill=manual.colors["Late"])+ 
      scale_y_continuous(name="# Genomic segments") + 
        scale_x_continuous(name="Mutations per Mb", limits = c(-0.05*(max(to.plot$Mean)), max.mutation.time)) +
      ggtitle(i)
    
    print(p)
    
    p <- ggplot(to.plot[to.plot$Timing=="Pre-CNV",], aes(x=Mean, fill=Segment)) + geom_histogram( binwidth = binwidth)+ 
      scale_fill_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot[to.plot$Timing=="Pre-CNV",])))+
      geom_density(data=to.plot[to.plot$Timing=="Post-CNV",], linetype=2, aes(x=Mean), inherit.aes = F) +
      scale_y_continuous(name="# Genomic segments") + 
        scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$Mean)), max(to.plot$Mean)*1.1)) +
      ggtitle(i)
    
    print(p)
    
    to.plot. <- to.plot[to.plot$Timing=="Pre-CNV",]
    to.plot.$Segment <- factor(to.plot.$Segment, levels=unique(to.plot.$Segment))
    
    p <- ggplot(to.plot., aes(x=Mean, xmin=Min, xmax=Max, y=as.numeric(Segment)/nrow(to.plot.), col=Segment)) +
      scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot.))) + 
      geom_point() + geom_errorbarh(height=0)+ 
      geom_vline(data=data.frame(x=mutation.time.mrca[i,]$Mean), aes(xintercept=x/3.3/10^3, y=1), inherit.aes = F, col=manual.colors["Late"],
                 linetype=2) + 
      geom_ribbon(data=data.frame(xmin=rep(mutation.time.mrca[i,]$Min,2),
                                  xmax=rep(mutation.time.mrca[i,]$Max,2),
                                  y=c(0, nrow(to.plot[to.plot$Timing=="Pre-CNV",]))), 
                  aes(xmin=xmin/3.3/10^3, xmax=xmax/3.3/10^3,y=y),
                  inherit.aes = F, fill=manual.colors["Late"], col=NA, alpha=0.5) +
      geom_vline(data=data.frame(x=mutation.time.eca[i,]$Mean), aes(xintercept=x/3.3/10^3, y=1), inherit.aes = F, col=manual.colors["Early"],
                 linetype=2) + 
      geom_ribbon(data=data.frame(xmin=rep(mutation.time.eca[i,]$Min, 2),
                                  xmax=rep(mutation.time.eca[i,]$Max,2),
                                  y=c(0, nrow(to.plot[to.plot$Timing=="Pre-CNV",]))), 
                  aes( xmin=xmin/3.3/10^3, xmax=xmax/3.3/10^3, y=y),
                  height=0, inherit.aes = F, fill=manual.colors["Early"], col=NA, alpha=0.5) +
      scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$Mean)), max(to.plot$Max)*1.1)) +
      scale_y_continuous(limits=c(0, nrow(to.plot)))
    
    print(p)

    dev.off()
  }
  
  if(i == "NBE44"){
    pdf(paste0(panel.directory, "Figure_S3d.pdf"))
    
    p <- ggplot(to.plot[to.plot$Timing=="Post-CNV",], aes(x=Mean)) + geom_histogram( binwidth = binwidth, fill=manual.colors["Late"])+ 
      scale_y_continuous(name="# Genomic segments") + 
       scale_x_continuous(name="Mutations per Mb", limits = c(-0.05*(max(to.plot$Mean)), max.mutation.time)) +
      ggtitle(i)
    
    print(p)
    
    p <- ggplot(to.plot[to.plot$Timing=="Pre-CNV",], aes(x=Mean, fill=Segment)) + geom_histogram( binwidth = binwidth)+ 
      scale_fill_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot[to.plot$Timing=="Pre-CNV",])))+
      geom_density(data=to.plot[to.plot$Timing=="Post-CNV",], linetype=2, aes(x=Mean), inherit.aes = F) +
      scale_y_continuous(name="# Genomic segments") + 
        scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$Mean)), max(to.plot$Mean)*1.1)) +
      ggtitle(i)
    
    print(p)
    
    to.plot. <- to.plot[to.plot$Timing=="Pre-CNV",]
    to.plot.$Segment <- factor(to.plot.$Segment, levels=unique(to.plot.$Segment))
    
    p <- ggplot(to.plot., aes(x=Mean, xmin=Min, xmax=Max, y=as.numeric(Segment)/nrow(to.plot.), col=Segment)) +
      scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot.))) + 
      geom_point() + geom_errorbarh(height=0)+ 
      geom_vline(data=data.frame(x=mutation.time.mrca[i,]$Mean), aes(xintercept=x/3.3/10^3, y=1), inherit.aes = F, col=manual.colors["Late"],
                 linetype=2) + 
      geom_ribbon(data=data.frame(xmin=rep(mutation.time.mrca[i,]$Min,2),
                                  xmax=rep(mutation.time.mrca[i,]$Max,2),
                                  y=c(0, nrow(to.plot[to.plot$Timing=="Pre-CNV",]))), 
                  aes(xmin=xmin/3.3/10^3, xmax=xmax/3.3/10^3,y=y),
                  inherit.aes = F, fill=manual.colors["Late"], col=NA, alpha=0.5) +
           scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$Mean)), max(to.plot$Max)*1.1)) +
      scale_y_continuous(limits=c(0, nrow(to.plot)))
    
    print(p)   
    
    dev.off()
  }
  
  if(i == "NBE99"){
    pdf(paste0(panel.directory, "Figure_S3a.pdf"))
    
    p <- ggplot(to.plot[to.plot$Timing=="Post-CNV",], aes(x=Mean)) + geom_histogram( binwidth = binwidth, fill=manual.colors["Late"])+ 
      scale_y_continuous(name="# Genomic segments") + 
         scale_x_continuous(name="Mutations per Mb", limits = c(-0.05*(max(to.plot$Mean)), max.mutation.time)) +
      ggtitle(i)
    
    print(p)
    
    p <- ggplot(to.plot[to.plot$Timing=="Pre-CNV",], aes(x=Mean, fill=Segment)) + geom_histogram( binwidth = binwidth)+ 
      scale_fill_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot[to.plot$Timing=="Pre-CNV",])))+
      geom_density(data=to.plot[to.plot$Timing=="Post-CNV",], linetype=2, aes(x=Mean), inherit.aes = F) +
      scale_y_continuous(name="# Genomic segments") + 
      scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$Mean)), max(to.plot$Mean)*1.1)) +
      ggtitle(i)
    
    print(p)
    
    to.plot. <- to.plot[to.plot$Timing=="Pre-CNV",]
    to.plot.$Segment <- factor(to.plot.$Segment, levels=unique(to.plot.$Segment))
    
    p <- ggplot(to.plot., aes(x=Mean, xmin=Min, xmax=Max, y=as.numeric(Segment)/nrow(to.plot.), col=Segment)) +
      scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot.))) + 
      geom_point() + geom_errorbarh(height=0)+ 
      geom_vline(data=data.frame(x=mutation.time.mrca[i,]$Mean), aes(xintercept=x/3.3/10^3, y=1), inherit.aes = F, col=manual.colors["Late"],
                 linetype=2) + 
      geom_ribbon(data=data.frame(xmin=rep(mutation.time.mrca[i,]$Min,2),
                                  xmax=rep(mutation.time.mrca[i,]$Max,2),
                                  y=c(0, nrow(to.plot[to.plot$Timing=="Pre-CNV",]))), 
                  aes(xmin=xmin/3.3/10^3, xmax=xmax/3.3/10^3,y=y),
                  inherit.aes = F, fill=manual.colors["Late"], col=NA, alpha=0.5) +
      scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$Mean)), max(to.plot$Max)*1.1)) +
      scale_y_continuous(limits=c(0, nrow(to.plot)))
    
    print(p)    

    dev.off()
  }
  
}


##############################################################################################################################################
##### Figure S2b: Compare ECA and MRCA between primary and relapse tumor

mutation.time.eca[earliest.mutation.time$Sample,]$Mean <- earliest.mutation.time$Mean
mutation.time.eca[earliest.mutation.time$Sample,]$Min <- earliest.mutation.time$Min
mutation.time.eca[earliest.mutation.time$Sample,]$Max <- earliest.mutation.time$Max

colnames(sample.information.discovery)[which(colnames(sample.information.discovery)=="ECA")] <- "ECA.exists"
to.plot <- cbind(sample.information.discovery, data.frame(MRCA=mutation.time.mrca[rownames(sample.information.discovery),]$Mean/3.3/10^3,
                                                    ECA=mutation.time.eca[rownames(sample.information.discovery),]$Mean/3.3/10^3,
                                                    MRCA.upper=mutation.time.mrca[rownames(sample.information.discovery),]$Max/3.3/10^3,
                                                    MRCA.lower=mutation.time.mrca[rownames(sample.information.discovery),]$Min/3.3/10^3,
                                                    ECA.upper=mutation.time.eca[rownames(sample.information.discovery),]$Max/3.3/10^3,
                                                    ECA.lower=mutation.time.eca[rownames(sample.information.discovery),]$Min/3.3/10^3))

to.plot$Time <- sample.information.discovery$Sample.type
to.plot$Time[to.plot$Time=="Relapse 3"] <- "Relapse"
to.plot$Time[to.plot$Time %in% c("Relapse tumor", "Relapse metastasis")] <- "Relapse"

to.plot$Time <- factor(to.plot$Time, levels=c("Primary", "Metastasis", "Relapse"))


addWorksheet(wb.s2, "2b_ECA")
writeData(wb.s2, "2b_ECA", to.plot[,c("ECA", "ECA.lower", "ECA.upper", "Rounded.ploidy", "Time")])
addWorksheet(wb.s2, "2b_MRCA")
writeData(wb.s2, "2b_MRCA", to.plot[,c("MRCA", "MRCA.lower", "MRCA.upper", "Rounded.ploidy", "Time")])


pdf(paste0(panel.directory, "Figure_S2b.pdf"), width=4, height=4, useDingbats = F)

p <- ggplot(to.plot[to.plot$ECA!=to.plot$MRCA,], aes(x=Time, y=ECA, col=as.character(Rounded.ploidy)))+geom_quasirandom() +
  scale_color_manual(values=c("2" = "black", "3" = "orange", "4" = "firebrick")) +
  scale_x_discrete( breaks = c("Primary FALSE", "Primary TRUE", "Metastasis FALSE"),
                    labels=c("Primary untreated", "Primary treated", "Metastasis"))+ 
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1,
         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_y_continuous(name = "SSNVs/Mb at ECA")

print(p)


chars <- capture.output(wilcox.test(to.plot[to.plot$Time=="Primary" & to.plot$ECA!= to.plot$MRCA,]$ECA, 
                                    to.plot[to.plot$Time=="Metastasis" & to.plot$ECA!= to.plot$MRCA,]$ECA, conf.int = T))

writeData(wb.s2, sheet = "2b_ECA", chars, startCol = 8)


chars <- capture.output(wilcox.test(to.plot[to.plot$Time=="Primary" & to.plot$ECA!= to.plot$MRCA,]$ECA, 
                                    to.plot[to.plot$Time=="Relapse" & to.plot$ECA!= to.plot$MRCA,]$ECA, conf.int = T))

writeData(wb.s2, sheet = "2b_ECA", chars, startCol = 10)

chars <- capture.output(wilcox.test(to.plot[to.plot$Time=="Relapse" & to.plot$ECA!= to.plot$MRCA,]$ECA,
                                    to.plot[to.plot$Time=="Metastasis" & to.plot$ECA!= to.plot$MRCA,]$ECA, conf.int = T))

writeData(wb.s2, sheet = "2b_ECA", chars, startCol = 12)


ggplot(to.plot, aes(x=Time, y=MRCA, col=as.character(Rounded.ploidy)))+geom_quasirandom() +
  scale_color_manual(values=c("2" = "black", "3" = "orange", "4" = "firebrick")) +
  scale_x_discrete( breaks = c("Primary FALSE", "Primary TRUE", "Metastasis FALSE"),
                    labels=c("Primary untreated", "Primary treated", "Metastasis"))+ 
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio=1,
         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_y_continuous(name = "SSNVs/Mb at MRCA")


chars <- capture.output(wilcox.test(to.plot[to.plot$Time=="Primary",]$MRCA, to.plot[to.plot$Time=="Metastasis",]$MRCA, conf.int = T))

writeData(wb.s2, sheet = "2b_MRCA", chars, startCol = 8)

chars <- capture.output(wilcox.test(to.plot[to.plot$Time=="Primary",]$MRCA, to.plot[to.plot$Time=="Relapse",]$MRCA, conf.int = T))

writeData(wb.s2, sheet = "2b_MRCA", chars, startCol = 10)

chars <- capture.output(wilcox.test(to.plot[to.plot$Time=="Metastasis",]$MRCA, to.plot[to.plot$Time=="Relapse",]$MRCA, conf.int = T))

writeData(wb.s2, sheet = "2b_MRCA", chars, startCol = 12)

dev.off()


##############################################################################################################################################
##### Figure S3b: Early MRCA with and without ECA; discovery cohort
cutpoint <- 0.05
pearly <- list()

max.mutation.time.primary <- max(mutation.time.mrca[primary.tumors.discovery,]$Max)

for(ECA.exists in c(T, F)){
  
  subset=sample.information.discovery[ sample.information.discovery$Sample.type %in% c("Primary", "Metastasis") &
                                   sample.information.discovery$ECA.exists==ECA.exists &
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
  
  
  
  if(ECA.exists){
    panel="S3b_ECA"
  }else{
    panel="S3b_no_ECA"
  }
  
  addWorksheet(wb.s3, panel)
  writeData(wb.s3, panel, to.plot[,c("ECA", "ECA.lower", "ECA.upper", "MRCA", "MRCA.lower", "MRCA.upper")])
  
  
  
  pearly[[length(pearly)+1]] <- ggplot(data = to.plot[order(to.plot$MRCA),],
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
    pearly[[length(pearly)]] <- pearly[[length(pearly)]] + 
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
  
}


pdf(paste0(panel.directory, "Figure_S3b.pdf"), width=9, height=9, useDingbats = F)

figure <- ggarrange(plotlist=pearly, nrow=3, ncol=3)
annotate_figure(figure, top="Primary tumor / Metastasis")

dev.off()


##########################################################################################################################################
saveWorkbook(wb, file = paste0(panel.directory, "Source_data_Fig.2.xlsx"), overwrite=T)

saveWorkbook(wb.s2, file = paste0(panel.directory, "Source_data_Fig.S2.xlsx"), overwrite=T)
saveWorkbook(wb.s3, file = paste0(panel.directory, "Source_data_Fig.S3.xlsx"), overwrite=T)


##############################################################################################################################################
## All tumors: density plots per genomic segment. Like Fig. 2d-i but for all tumors of the discovery set.

p.all.primary.tumors <- list()
p.all.relapse.tumors <- list()
#p.mrca.densities.all.tumors <- list()

## Set the range of the samples of interest
#max.mutation.time <- 0.7
max.mutation.time.primary.tumors <- max(unlist(lapply(mutation.time[primary.tumors.discovery], function(x){x$Mean}))/3.3/10^3)
max.mutation.time.relapse.tumors <- max(unlist(lapply(mutation.time[relapse.tumors.discovery], function(x){x$Mean}))/3.3/10^3)

## Moreover quantify the mutation density of early and late clonal mutations
mutational.density.post <- c()
mutational.density.pre <- c()

source.data <- c()

for(i in tumors.discovery){

  binwidth = (max(mutation.time[[i]]$Mean/3.3/10^3) - min(mutation.time[[i]]$Mean/3.3/10^3))/20
  
  to.plot <- mutation.time[[i]]
  to.plot$Mean <- to.plot$Mean/3.3/10^3
  to.plot$Min <- to.plot$Min/3.3/10^3
  to.plot$Max <- to.plot$Max/3.3/10^3
  
  to.plot$Timing <- sapply(as.character(to.plot$Segment), function(x){
                          x <- strsplit(x, split="_")[[1]]
                          if(is.na(x[4])){
                            return("Post-CNV")
                          }
                          if(x[4]=="I"){
                            return("Post-CNV")
                          }else{
                            return("Pre-CNV")
                          }
                        })
  
  to.plot$Segment <- factor(to.plot$Segment, levels=unique(to.plot$Segment))
  to.plot$Tumor <- i
  
  source.data <- rbind(source.data, to.plot)
  
  
  
  if(i %in% primary.tumors.discovery){
    p.all.primary.tumors[[length(p.all.primary.tumors)+1]] <-  ggplot(to.plot, aes(x=Mean, fill=Timing)) + geom_histogram( binwidth = 0.05, alpha=0.75)+ 
      scale_fill_manual(values=c("Pre-CNV"=unname(manual.colors["Early"]), "Post-CNV" = unname(manual.colors["Late"])))+
      scale_y_continuous(name="# Genomic segments") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      scale_x_continuous(name="Mutations per Mb", limits=c(-0.075, max.mutation.time.primary.tumors)) +
      ggtitle(sample.information.discovery[i,]$Evolution_paper_Id)
  }
  
  if(i %in% relapse.tumors.discovery){
    p.all.relapse.tumors[[length(p.all.relapse.tumors)+1]] <-  ggplot(to.plot, aes(x=Mean, fill=Timing)) + geom_histogram( binwidth = 0.05, alpha=0.75)+ 
      scale_fill_manual(values=c("Pre-CNV"=unname(manual.colors["Early"]), "Post-CNV" = unname(manual.colors["Late"])))+
      scale_y_continuous(name="# Genomic segments") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
      scale_x_continuous(name="Mutations per Mb", limits=c(-0.075, max.mutation.time.relapse.tumors)) +
      ggtitle(sample.information.discovery[i,]$Evolution_paper_Id)
  }
  
  
  
  
  pre.clonal.mutations <- which(sapply(as.character(mutation.time[[i]]$Segment), function(x){
    x <- strsplit(x, split="_")[[1]]
    if(is.na(x[4])){
      return(0)
    }
    if(x[4]=="I"){
      return(0)
    }else{
      return(1)
    }
  })==1)
  
  post.clonal.mutations <- which(sapply(as.character(mutation.time[[i]]$Segment), function(x){
    x <- strsplit(x, split="_")[[1]]
    if(is.na(x[4])){
      return(0)
    }
    if(x[4]=="I"){
      return(0)
    }else{
      return(1)
    }
  })==0)
  
  mutational.density.pre <- c(mutational.density.pre, median(mutation.time[[i]][pre.clonal.mutations,]$Mean)/3.3/10^3)
  mutational.density.post <- c(mutational.density.post, median(mutation.time[[i]][post.clonal.mutations,]$Max)/3.3/10^3)
  names(mutational.density.post)[length(mutational.density.post)] <- i
  names(mutational.density.pre)[length(mutational.density.post)] <- i
  
}
names(p.all.primary.tumors) <- sample.information.discovery[primary.tumors.discovery, ]$Evolution_paper_Id
names(p.all.relapse.tumors) <- sample.information.discovery[relapse.tumors.discovery, ]$Evolution_paper_Id


pdf(paste0(panel.directory, "Mutation_densities_all_tumors.pdf"), width=10, height=15)

p.primary.tumors.triploid.eca <- sample.information.discovery[sample.information.discovery$ECA.exists==T & sample.information.discovery$Rounded.ploidy==3 & 
                                                          sample.information.discovery$Sample.type %in% c("Primary", "Metastasis"),]$Evolution_paper_Id
p.primary.tumors.di.tetra.eca<- sample.information.discovery[sample.information.discovery$ECA.exists==T & sample.information.discovery$Rounded.ploidy!=3 & 
                                                         sample.information.discovery$Sample.type %in% c("Primary", "Metastasis"),]$Evolution_paper_Id
p.primary.tumors.triploid.no.eca<- sample.information.discovery[sample.information.discovery$ECA.exists==F & sample.information.discovery$Rounded.ploidy==3 & 
                                                            sample.information.discovery$Sample.type %in% c("Primary", "Metastasis"),]$Evolution_paper_Id
p.primary.tumors.di.tetra.no.eca<- sample.information.discovery[sample.information.discovery$ECA.exists==F & sample.information.discovery$Rounded.ploidy!=3 & 
                                                            sample.information.discovery$Sample.type %in% c("Primary", "Metastasis"),]$Evolution_paper_Id

ggarrange(plotlist=p.all.primary.tumors[p.primary.tumors.triploid.eca], nrow=8, ncol=4, common.legend = T)
ggarrange(plotlist=p.all.primary.tumors[p.primary.tumors.triploid.no.eca], nrow=8, ncol=4, common.legend = T)
ggarrange(plotlist=p.all.primary.tumors[p.primary.tumors.di.tetra.eca], nrow=8, ncol=4, common.legend = T)
ggarrange(plotlist=p.all.primary.tumors[p.primary.tumors.di.tetra.no.eca], nrow=8, ncol=4, common.legend = T)

ggarrange(plotlist=p.all.relapse.tumors, nrow=8, ncol=4, common.legend = T)

dev.off()



#########################################################################################################################################
###### all tumors

## manually set the range 
p.all.tumors.post.cnv <- list()
p.all.tumors.pre.cnv <- list()
p.all.timeline <- list()

## workbook for source data
timing.info <- c()
for(i in names(mutation.time)){

  binwidth = (max(mutation.time[[i]]$Mean/3.3/10^3) - min(mutation.time[[i]]$Mean/3.3/10^3))/20
  
  to.plot <- mutation.time[[i]]
  to.plot$Mean <- to.plot$Mean/3.3/10^3
  to.plot$Min <- to.plot$Min/3.3/10^3
  to.plot$Max <- to.plot$Max/3.3/10^3
  
  to.plot$Timing <- sapply(as.character(to.plot$Segment), function(x){
    x <- strsplit(x, split="_")[[1]]
    if(is.na(x[4])){
      return("Post-CNV")
    }
    if(x[4]=="I"){
      return("Post-CNV")
    }else{
      return("Pre-CNV")
    }
  })
  
 
  to.plot$Segment <- factor(to.plot$Segment, levels=unique(to.plot$Segment))
 
    
    p.all.tumors.post.cnv[[i]] <- eval(substitute(ggplot(to.plot[to.plot$Timing=="Post-CNV",], aes(x=Mean)) + 
      geom_histogram( binwidth = binwidth, fill=manual.colors["Late"])+ 
      scale_y_continuous(name="# Genomic segments") + 
      scale_x_continuous(name="Mutations per Mb", limits = c(-0.05*(max(to.plot$Mean)), max(to.plot$Mean)*1.1)) +
      ggtitle(i), list(i=i, to.plot=to.plot)))
    

    p.all.tumors.pre.cnv[[i]] <- eval(substitute(ggplot(to.plot[to.plot$Timing=="Pre-CNV",], aes(x=Mean)) + 
      geom_histogram( binwidth = binwidth, fill=manual.colors["Early"])+ 
      # scale_fill_manual(values=brewer.pal(9, "PRGn")[c(9)])+
      #scale_fill_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot[to.plot$Timing=="Pre-CNV",])))+
      geom_density(data=to.plot[to.plot$Timing=="Post-CNV",], linetype=2, aes(x=Mean), inherit.aes = F) +
      scale_y_continuous(name="# Genomic segments") + 
      scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$Mean)), max(to.plot$Mean)*1.1)) +
      ggtitle(i), list(i=i, to.plot=to.plot)))
    

    to.plot. <- to.plot[to.plot$Timing=="Pre-CNV",]
    to.plot.$Segment <- factor(to.plot.$Segment, levels=unique(to.plot.$Segment))
    
    ## add significance info
    
    to.plot.$Time <- sapply(to.plot.$Segment, function(x){
      x <- as.character(x)
      if(x %in% mrca.eca[[i]]$gains.at.mrca & x %in% mrca.eca[[i]]$gains.at.mrca.conforming.eca){
        "ECA/MRCA"
      }else if(x %in% mrca.eca[[i]]$gains.at.mrca){
        "MRCA"
      }else if(x %in% mrca.eca[[i]]$gains.uniquely.mapped.to.eca){
        "ECA"
      }else if(x %in% mrca.eca[[i]]$gains.at.earliest.time){
        "< ECA"
      }else if(x %in% mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca/3.3/10^3 & 
               to.plot.[to.plot.$Segment==x,]$Min > mrca.eca[[i]]$mutation.time.mrca/3.3/10^3 & 
               to.plot.[to.plot.$Segment==x,]$Min <= mrca.eca[[i]]$mutation.time.mrca.upper/3.3/10^3){
        "MRCA"
      }else if(x %in% mrca.eca[[i]]$gains.not.maping.to.eca.or.mrca & !is.na(mrca.eca[[i]]$mutation.time.eca) &&
                to.plot.[to.plot.$Segment==x,]$Mean > mrca.eca[[i]]$mutation.time.eca/3.3/10^3 && 
                to.plot.[to.plot.$Segment==x,]$Mean < mrca.eca[[i]]$mutation.time.mrca/3.3/10^3){
        "> ECA, < MRCA"
      }else{
        "> MRCA"
      }
    })
    
    if(nrow(to.plot.)>0){
      timing.info <- rbind(timing.info,
                           cbind(to.plot., data.frame(Sample=i)))
    }

    
    to.plot.$Segment <- sapply(to.plot.$Segment, function(x){
      x <- as.character(x)
      x <- strsplit(x, split="_")[[1]]
      if(x[3]=="disomic"){
        x[3] <- "2N"
      }else if(x[3]=="trisomic"){
        x[3] <- "3N"
      }else if(x[3]=="monosomic"){
        x[3] <- "1N"
      }else{
        x[3] <- "4N"
      }
      paste(x[c(2,3,4)], collapse="-")
    })
    
    to.plot.$Segment <- factor(to.plot.$Segment, levels=unique(as.character(to.plot.$Segment)))
    
    
    ## if there is only one segment, place it in the center, else distribute equally
    p.all.timeline[[i]] <- eval(substitute(ggplot(to.plot., aes(x=Mean, xmin=Min, xmax=Max, y=as.numeric(Segment) - 
                                                  0.5, col = Segment)) +
    # scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(nrow(to.plot.))) + 
      geom_point(aes(shape=Time)) + geom_errorbarh(height=0)+ 
     # geom_text(aes(label=str_wrap(Segment,12), y = as.numeric(Segment) - 0.25), size=6) +
      geom_vline(data=data.frame(x=mutation.time.mrca[i,]$Mean), aes(xintercept=x/3.3/10^3, y=1), inherit.aes = F, col=manual.colors["Late"],
                 linetype=2) + 
      geom_ribbon(data=data.frame(xmin=rep(mutation.time.mrca[i,]$Min,2),
                                  xmax=rep(mutation.time.mrca[i,]$Max,2),
                                  y=c(0, max(1,nrow(to.plot.)))), 
                  aes(xmin=xmin/3.3/10^3, xmax=xmax/3.3/10^3,y=y),
                  inherit.aes = F, fill=manual.colors["Late"], col=NA, alpha=0.5) +
      geom_vline(data=data.frame(x=mutation.time.eca[i,]$Mean), aes(xintercept=x/3.3/10^3, y=1), inherit.aes = F, col=manual.colors["Early"],
                 linetype=2) + 
      geom_ribbon(data=data.frame(xmin=rep(mutation.time.eca[i,]$Min, 2),
                                  xmax=rep(mutation.time.eca[i,]$Max,2),
                                  y=c(0, nrow(to.plot.))), 
                  aes( xmin=xmin/3.3/10^3, xmax=xmax/3.3/10^3, y=y),
                  height=0, inherit.aes = F, fill=manual.colors["Early"], col=NA, alpha=0.5) +
      scale_x_continuous(name="Mutations per Mb", limits =  c(-0.05*(max(to.plot$Mean)), max(to.plot$Max)*1.1)) +
      scale_y_continuous(limits=c(0, max(1,nrow(to.plot.)))) + theme(axis.line.y=element_blank(), axis.text.y=element_blank(),
                                                 axis.ticks.y=element_blank(), axis.title.y=element_blank()) +
      theme(legend.text = element_text(size=6)),
    list(to.plot=to.plot, i=i, to.plot.=to.plot.)))
  
}

save(p.all.tumors.post.cnv, file=paste0(rdata.directory, "Post_CNV_dens_distr_all_tumors.RData"))

save(p.all.tumors.pre.cnv, file=paste0(rdata.directory, "Pre_CNV_dens_distr_all_tumors.RData"))

save(p.all.timeline, file=paste0(rdata.directory, "Timeline_all_tumors.RData"))