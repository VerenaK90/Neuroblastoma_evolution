##############################################################################################################################################
## Plot the results per sample
##############################################################################################################################################
source("Settings.R")
load(paste0(rdata.directory, "Purity_ploidy.RData"))
detach("package:ggbio", unload = TRUE)
library(cowplot)
library(wesanderson)

load(paste0(rdata.directory, "Purity_ploidy.RData"))
load(paste0(rdata.directory, "VAFs_all_tumors.RData"))
load(paste0(rdata.directory, "Pre_CNV_dens_distr_all_tumors.RData"))
load(paste0(rdata.directory, "Post_CNV_dens_distr_all_tumors.RData"))
load(paste0(rdata.directory, "Timeline_all_tumors.RData"))


## all mutations
sigs.discovery <- read.delim(paste0(signature.directory, "/discovery/MMSig_all_discovery.tsv"))
## sum up SBS1, SBS5 and SBS40 to a "clock-like" signature
sigs.discovery$Clock <- sigs.discovery$SBS1 + sigs.discovery$SBS5 + sigs.discovery$SBS40
sigs.discovery$Sample <- rownames(sigs.discovery)

## clonal SNVs only
clonal.sigs.discovery <- read.delim(paste0(signature.directory, "/discovery/MMSig_output_discovery_C.tsv"))
clonal.sigs.discovery$Clock <- clonal.sigs.discovery$SBS1 + clonal.sigs.discovery$SBS5 + clonal.sigs.discovery$SBS40
clonal.sigs.discovery$Sample <- rownames(clonal.sigs.discovery)

## subclonal SNVs only
subclonal.sigs.discovery <- read.delim(paste0(signature.directory, "/discovery/MMSig_output_discovery_SC.tsv"))
subclonal.sigs.discovery$Clock <- subclonal.sigs.discovery$SBS1 + subclonal.sigs.discovery$SBS5 + subclonal.sigs.discovery$SBS40
subclonal.sigs.discovery$Sample <- rownames(subclonal.sigs.discovery)

## validation cohort
sigs.validation <- read.delim(paste0(signature.directory, "/validation/MMSig_all_validation.tsv"))
## sum up SBS1, SBS5 and SBS40 to a "clock-like" signature
sigs.validation$Clock <- sigs.validation$SBS1 + sigs.validation$SBS5 + sigs.validation$SBS40
sigs.validation$Sample <- rownames(sigs.validation)


load(paste0(rdata.directory, "Sig_colors.RData"))
##############################################################################################################################################

summary.all.tumors <- list()

for(i in c(tumors.discovery, tumors.validation)){

  if(i %in% tumors.discovery){
    data.directory <- data.directory.discovery
  }else{
    data.directory <- data.directory.validation
  }
  
  # add a common title
  
  id <- ifelse(i %in% tumors.discovery, sample.information.discovery[i,]$Evolution_paper_Id,
               sample.information.validation[i,]$Evolution_paper_ID)
  
  title <- ggdraw() + 
    draw_label(
      id,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 2)
    )
  
  ## First and second row: quality measures (karyogram, ploidy and copy number profile), as well as mutational signatures

  load(paste0(data.directory, i, "/CNAqc/Quality_check.RData"))
  ## ploidy estimate compariston between DNA-index and ACEseq, if information is available
  ploidy.aceseq <- ploidy[i]
  if(i %in% tumors.discovery & !is.na(as.numeric(sample.information.discovery[i,]$`DNA-INDEX`))){
    ploidy.dna.index <- 2*as.numeric(sample.information.discovery[i,]$`DNA-INDEX`)
    
    to.plot <- data.frame(Method=c("Aceseq", "DNA-Index"), Ploidy=c(ploidy.aceseq, ploidy.dna.index))
    
    ploidy.plot <- ggplot(to.plot, aes(x=Method, y=Ploidy)) + geom_col() + scale_y_continuous(name="Average ploidy estimate")
  }else{
    to.plot <- data.frame(Method=c("Aceseq"), Ploidy=ploidy.aceseq)
    
    ploidy.plot <- ggplot(to.plot, aes(x=Method, y=Ploidy)) + geom_col() + scale_y_continuous(name="Average ploidy estimate")
  }
  
  panel.count.first.row <- 2
  

  first_row <- plot_grid(segments, ploidy.plot, ncol=2, labels=c("a", "b"), label_size = 12,
                       rel_widths = c(6, 1), label_fontfamily = "Helvetica", label_colour = "black", align = "hv")
  

  ##mutational signatures; all, clonal and subclonal
  
  if(i %in% tumors.discovery){
    to.plot <- rbind(cbind(melt(sigs.discovery[sigs.discovery$Sample==i,-which(colnames(sigs.discovery)=="mutations"),drop=F], id.vars = "Sample",
                                variable.name = "Signature", value.name = "Fraction of SNVs"), data.frame(Type="All")),
                     cbind(melt(clonal.sigs.discovery[clonal.sigs.discovery$Sample==i,-which(colnames(clonal.sigs.discovery)=="mutations"),drop=F], id.vars = "Sample",
                                variable.name = "Signature", value.name = "Fraction of SNVs"), data.frame(Type="Clonal")),
                     cbind(melt(subclonal.sigs.discovery[subclonal.sigs.discovery$Sample==i,-which(colnames(subclonal.sigs.discovery)=="mutations"),drop=F], id.vars = "Sample",
                                variable.name = "Signature", value.name = "Fraction of SNVs"), data.frame(Type="Subclonal")))
    
    to.plot <- to.plot[!to.plot$Signature %in% c("SBS1", "SBS5", "SBS40"),]
    p1 <- ggplot(to.plot, aes(x=Type, y=`Fraction of SNVs`, fill=Signature)) + geom_col() + scale_fill_manual(values=sig.colors) +
      scale_x_discrete(name="") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      theme(legend.key.size = unit(0.5, "cm"))
  }else{
    to.plot <- cbind(melt(sigs.validation[sigs.validation$Sample==i,-which(colnames(sigs.validation)=="mutations"),drop=F], id.vars = "Sample",
                                variable.name = "Signature", value.name = "Fraction of SNVs"), data.frame(Type="All"))
    
    to.plot <- to.plot[!to.plot$Signature %in% c("SBS1", "SBS5", "SBS40"),]
    p1 <- ggplot(to.plot, aes(x=Type, y=`Fraction of SNVs`, fill=Signature)) + geom_col() + scale_fill_manual(values=sig.colors) +
      scale_x_discrete(name="") 
  }
  
  panel.count.second.row <- panel.count.first.row + 2
  
  if(length(q.1)==0){
    q.1 <- ggplot() + theme(axis.line=element_blank())
  }
  
  second_row <- plot_grid(plot_grid(plotlist=q.1, nrow=1, ncol=4),
                          plot_grid(p1, align = "hv"),
                          ncol=2, rel_widths = c(3, 1), rel_heights = c(1, 0.1),
                          labels = letters[(panel.count.first.row+1) : panel.count.second.row],
                          label_size = 12, label_fontfamily = "Helvetica", label_colour = "black")
  
  
  ## Third row: mutation density at ECA, MRCA and timeline
  
  panel.count.third.row <- panel.count.second.row + 1
  
  third_row <- plot_grid(p.all.tumors.post.cnv[[i]] + ggtitle(""), p.all.tumors.pre.cnv[[i]] + ggtitle(""), p.all.timeline[[i]] + ggtitle("") + 
                           theme(legend.key.size = unit(0.25, "cm")), 
                          labels=letters[(panel.count.second.row + 1):panel.count.third.row], ncol=3, label_size = 12, rel_widths = c(1,1,3),
                          label_fontfamily = "Helvetica", label_colour = "black", align = "hv")
  
  
  ## Fourth row: Mobster and model fits
  
  p.evol <- list()
  
  if(file.exists(paste0(data.directory.discovery, "/", i, "/Mobster_fits.RData"))){
    load(paste0(data.directory.discovery, "/", i, "/Mobster_fits.RData"))
    
    p.evol[[1]] <- plot(fit$best) + ggtitle("") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=8, color="black"),
                                                        panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.75),
                                                        axis.ticks = element_line(color = "black", size=0.75),
                                                        axis.text = element_text(size=8, color="black"))
    
    panel.count.fourth.row <- panel.count.third.row + 1
    
  }else{
    panel.count.fourth.row <- panel.count.third.row
    p.evol[[1]] <- ggplot() + theme_transparent()       
  }
  
  ##model fits 
  
  if(i %in% rownames(sample.information.discovery[sample.information.discovery$Use.to.fit.tumor.expansion=="T" ,])){
    load(paste0(data.directory.discovery, i, "/pyABC/Plot_fit.RData"))
    
    for(j in 1:length(p)){
      p.evol[[j+1]] <- p[[j]]
    }
    
    fits <- read.csv(paste0(data.directory.discovery, i, "/pyABC/", i, ".csv"))
    
    p.evol[[length(p.evol)+1]] <- ggplot(data.frame(Median= median(fits$par_mu/(1-fits$par_delta)),
                                                    Max = hdi(fits$par_mu/(1-fits$par_delta), credMass = 0.8)[1],
                                                    Min = hdi(fits$par_mu/(1-fits$par_delta), credMass = 0.8)[2]), 
                                         aes(x="", y=Median, ymin=Min, ymax=Max)) + geom_pointrange() +
      scale_y_continuous(limits=c(0, hdi(fits$par_mu/(1-fits$par_delta), credMass = 0.8)[2]*1.1),
                         name="# SNVs per effective division") +
      scale_x_discrete(name="Mutation rate")
    
    panel.count.fourth.row <- panel.count.fourth.row + 2
    
    panel.labels.fourth.row <- c(letters[(panel.count.third.row+1):(panel.count.third.row+2)], 
                               rep("", length(p.evol) - panel.count.fourth.row + panel.count.third.row),
                               letters[(panel.count.third.row+3)])
    

    
  }else{
    if(panel.count.third.row==panel.count.fourth.row){
      panel.labels.fourth.row <- ""
    }else{
      panel.labels.fourth.row <- c(letters[panel.count.fourth.row], rep("", 5 - panel.count.fourth.row + panel.count.third.row))
      p.evol[[2]] <- ggplot() + theme_transparent()           
      p.evol[[3]] <- ggplot() + theme_transparent()                                  
      p.evol[[4]] <- ggplot()  + theme_transparent()                             
      p.evol[[5]] <- ggplot() + theme_transparent()                            
      
    }
  }
  
  fourth_row <- plot_grid(plotlist=p.evol, 
                         labels= panel.labels.fourth.row,
                         ncol=length(p.evol), label_size = 12, 
                         label_fontfamily = "Helvetica", label_colour = "black", align = "hv")
  
  summary.all.tumors[[i]] <- plot_grid(title, first_row, second_row, third_row, fourth_row, 
                                       ncol = 1, nrow=5, rel_heights = c(0.1, 0.75, 0.75, 1, 1), align = "hv")
    
}

pdf(paste0(output.directory, "Summary_all_patients.pdf"), width=8.5, height = 11)

for(i in 1:length(summary.all.tumors)){
  print(summary.all.tumors[[i]])
}
dev.off()
