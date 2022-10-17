## Reproduce Fig. 1 - discovery data set only
##############################################################################################################################################
## load settings and libraries

source("Settings.R")
source(paste0(custom.script.directory, "Gains_and_losses.R"))
source(paste0(custom.script.directory, "Driver_mutations.R"))
source(paste0(custom.script.directory, "Oncoprint_Fig1.R"))
#source(paste0(custom.script.directory, "Mutational_signatures.R")) ## source or re-run

## store source data
## Figure 1
wb <- createWorkbook()
## Figure S1
wb.s <- createWorkbook()

## store figure panels
panel.directory <- paste0(output.directory, "Figure1/")

if(!dir.exists(panel.directory)){
  dir.create(panel.directory)
}

##############################################################################################################################################
## Figure 1a: overview of patient cohort

to.plot <- sample.information.discovery[,c("Telomere.maintenance.mechanism", "Sample.type", "Age", "ManualScore", "Stage")]
to.plot$ID <- rownames(to.plot)
## manipulate data slightly:
to.plot$Age.binary <- "<18 months"
to.plot$Age.binary[to.plot$Age >= 18*30]  <- ">=18 months"
to.plot$Age <- to.plot$Age.binary
to.plot$Age <- factor(to.plot$Age, levels=c("<18 months", ">=18 months"))

to.plot$Stage[to.plot$Stage==5] <- "4S"
to.plot$Stage <- factor(to.plot$Stage, levels=c("1", "4S", "2", "3", "4"))

##

to.plot$Sample.type[to.plot$Sample.type %in% c("Primary", "Metastasis")] <- "Primary"
to.plot$Sample.type <- factor(to.plot$Sample.type, levels=c("Primary", "Relapse"))
to.plot$Telomere.maintenance.mechanism <- factor(to.plot$Telomere.maintenance.mechanism, levels=c("MNA", "TERT", "ALT", "Multiple", "None"))
to.plot$Sample_pair <- "No pair"
to.plot$Sample_pair[to.plot$ID %in% c("NBE51", "NBE78")] <- "NBE51/NBE78"
to.plot$Sample_pair[to.plot$ID %in% c("NBE11", "NBE66")] <- "NBE11/NBE66"
to.plot$Sample_pair <- factor(to.plot$Sample_pair)


to.plot <- to.plot[order(to.plot$Telomere.maintenance.mechanism, decreasing = F),]
to.plot <- to.plot[order(to.plot$Age, decreasing = F),]
to.plot <- to.plot[order(to.plot$Sample.type, decreasing = F),]
to.plot <- to.plot[order(to.plot$Stage, decreasing = F),]

to.store <- to.plot
to.store$ID <- sample.information.discovery[to.plot$ID,]$Evolution_paper_Id

addWorksheet(wb, "a")
writeData(wb, "a", to.store)

pdf(paste0(panel.directory, "/Figure_1a.pdf"))

##separate heatmaps for high-risk, intermediate-risk, low-risk and the 2 sample pairs

ht_list <-  Heatmap(to.plot$Stage, name="Stage", col = stage.colors, width=unit(0.5, "cm"), 
                    split=to.plot[,c("Sample_pair", "ManualScore")], rect_gp = gpar(col="white", lwd=2)) +
  Heatmap(to.plot$Sample.type, name="Sample", col = time.colors, width=unit(0.5, "cm"),  rect_gp = gpar(col="white", lwd=2)) +
  Heatmap(as.numeric(to.plot$Age), name="Age", col = c("lightgrey", "black"),
          heatmap_legend_param = list(labels=levels(to.plot$Age), at=c(1,2)), width=unit(0.5, "cm"), rect_gp = gpar(col="white", lwd=2)) +
  Heatmap(to.plot$Telomere.maintenance.mechanism, name="TMM", col = telomere.colors, width=unit(0.5, "cm"),  rect_gp = gpar(col="white", lwd=2)) 

draw(ht_list, ht_gap=unit(0, "cm"))


dev.off()

## summary statistics:

## Sample type:
table(sample.information.discovery$Sample.type)

## 7 Metastasis, 60 Primary, 23 Relapse tumors, 10 Relapse metastasis

## Age at diagnosis:
median(sample.information.discovery$Age[sample.information.discovery$Stage %in% c(1, 2, 5)])/365
## 0.7y
median(sample.information.discovery$Age[sample.information.discovery$Stage ==3])/365
## 3.5y
median(sample.information.discovery$Age[sample.information.discovery$Stage ==4])/365
## 3.9y


##############################################################################################################################################
## Figure 1b: Purity/Ploidy overview

pdf(paste0(panel.directory, "/Figure_1b.pdf"), useDingbats = F, width=4, height=4)

addWorksheet(wb, "b")
writeData(wb, "b", sample.information.discovery[,c("Rounded.ploidy", "Purity")])

ggplot(sample.information.discovery, aes(x=Rounded.ploidy, y=Purity, group=Rounded.ploidy)) + geom_beeswarm() + scale_y_continuous(limits=c(0,1)) + 
  geom_boxplot(fill=NA)

dev.off()

median(sample.information.discovery$Purity)
## 0.88

##############################################################################################################################################
## Figure S1a: age distribution

pdf(paste0(panel.directory, "Figure_S1a.pdf"), width = 5, height=4, useDingbats = F)

to.plot <- sample.information.discovery[sample.information.discovery$Sample.type %in% c("Primary", "Metastasis"), c("Age", "Stage", "Telomere.maintenance.mechanism")] 
to.plot$Stage[to.plot$Stage=="5"] <- "4S"

addWorksheet(wb.s, "a")
writeData(wb.s, "a", to.plot) 

p <- ggplot(to.plot, aes(x=Age/365, col=Telomere.maintenance.mechanism)) + 
  stat_ecdf(size=2)+scale_color_manual(values=telomere.colors) + scale_y_continuous(name="Cumulative incidence") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(p)


p <- ggplot(to.plot, aes(x=Age/365, col=Stage)) + 
  stat_ecdf(size=2)+scale_color_manual(values=stage.colors) + scale_y_continuous(name="Cumulative incidence") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=10),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(p)

dev.off()


##############################################################################################################################################
## Figure S1b: Ploidy vs. DNA-index*2

pdf(paste0(panel.directory, "/Figure_S1b.pdf"), useDingbats = F, width=4, height=4)

load(paste0(rdata.directory, "Purity_ploidy.RData"))

to.plot <- data.frame(DNAIndex = as.numeric(sample.information.discovery[tumors.discovery,]$`DNA-INDEX`),
                      ACEseq = ploidy[tumors.discovery])

addWorksheet(wb.s, "b")
writeData(wb.s, "b", to.plot)

ggplot(to.plot, aes(x=ACEseq, y=DNAIndex*2)) + geom_point() + scale_y_continuous(limits=c(1,max(to.plot)*1.1)) +
  scale_x_continuous(limits=c(1, max(to.plot)*1.1)) + geom_abline(slope=1, intercept = 0, linetype=2)

cor.test(to.plot$DNAIndex*2, to.plot$ACEseq )
## 0.67

dev.off()


##############################################################################################################################################
## Figure S1c: Plot association between ploidy and stage

pdf(paste0(panel.directory, "Figure_S1c.pdf"))

to.plot <- data.frame(Stage = sample.information.discovery[tumors.discovery,]$Stage,
                      Ploidy = sample.information.discovery$Rounded.ploidy )
to.plot$Ploidy <- factor(to.plot$Ploidy, levels = c("2", "3", "4"))
to.plot$Stage <-  replace(to.plot$Stage, to.plot$Stage=="5", "4S")
to.plot$Stage <- factor(to.plot$Stage, levels = c("1", "4S", "2", "3", "4"))

addWorksheet(wb.s, "c")
writeData(wb.s, "c", to.plot)


to.plot <- to.plot %>% 
  group_by(Ploidy, Stage) %>% 
  summarise(count = n()) %>% 
  mutate(perc = count/sum(count))

ggplot(data = to.plot,
       aes(x=Ploidy, fill=Stage, y=perc*100)) + geom_bar(stat="identity", position="dodge", width=0.7) +
  scale_fill_manual(values =stage.colors) + scale_y_continuous(name="% Tumors") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

dev.off()

table(sample.information.discovery$Rounded.ploidy)
## 2: 55, 3: 33, 4: 12

##############################################################################################################################################
## Figure 1c: Plot numbers of chromosomes with gains and losses

## now plot the number of fragments that are clonally gained and the number of fragments that are subclonally gained

pdf(paste0(panel.directory, "Figure_1c.pdf"), useDingbats = F, width=2.5, height = 3)

## gains and losses, not sorted by subgroup

to.plot <- data.frame(Counts = c(chromosomes.with.clonal.gains, chromosomes.with.subclonal.gains,
                                 chromosomes.with.clonal.losses, chromosomes.with.subclonal.losses),
                      CNVType = c(rep("Gains", length(chromosomes.with.clonal.gains), length(chromosomes.with.subclonal.gains)*2),
                                  rep("Loss", length(chromosomes.with.clonal.losses), length(chromosomes.with.subclonal.losses)*2)), 
                      Clonality = c(rep(c("Clonal", "Subclonal"), each = length(chromosomes.with.clonal.gains)),
                                    rep(c("Clonal", "Subclonal"), each = length(chromosomes.with.clonal.losses))))

addWorksheet(wb, "c")
writeData(wb, "c", to.plot) 


ggplot(to.plot, aes(x=CNVType, y = Counts, fill = Clonality))+ geom_bar(stat="summary", fun.y="mean", position = "dodge", width=0.75)+
  geom_errorbar(stat="summary", position = "dodge", width=0.75) + scale_y_continuous(name="# Chromosomes")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=12),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c(Clonal = "black", Subclonal = "grey" )) + theme(legend.position = "bottom")

dev.off()

##############################################################################################################################################
## Figure S1d: Gains and losses across the cohort

pdf(paste0(panel.directory, "Figure_S1d.pdf"))

to.store <- rbind(data.frame(Chrom=seqnames(gain.granges), Start = start(gain.granges), End = end(gain.granges),
                       Number.of.cases = gain.granges$Coverage, Type = gain.granges$Type),
                  data.frame(Chrom=seqnames(loss.granges), Start = start(loss.granges), End = end(loss.granges),
                             Number.of.cases = loss.granges$Coverage, Type = loss.granges$Type))

addWorksheet(wb.s, "d")
writeData(wb.s, "d", to.store) 

plotGrandLinear(c(gain.granges, loss.granges), coord="genome", geom="bar", 
                ymax=c(gain.granges, loss.granges)$Coverage/length(tumors.discovery)*100, aes(y=Coverage/length(tumors.discovery)*100,
                                                                                    col=Type, fill=Type)) +  theme(panel.background = element_rect(fill = "white", color="black"), panel.grid.major = element_blank(),
                                                                                                                   panel.grid.minor = element_blank()) + geom_hline(yintercept = 0) +
  scale_color_manual(values = c("loss" = "#DA1F28", "gain"= "#B0D8A2")) +scale_fill_manual(values = c("loss" = "firebrick", "gain"= "darkgreen"))+
  scale_y_continuous("% Tumors")
dev.off()


##############################################################################################################################################
## Fig. 1d: Oncoprint


pdf(paste0(panel.directory, "/Figure_1d.pdf"), useDingbats = F, width=10, height = 7)

column_annotation <- data.frame(Ploidy = unlist(sapply(tumors.discovery, function(x){
  if(x %in% diploid.tumors.discovery){
    "2"
  }else if(x %in% triploid.tumors.discovery){
    "3"
  }else if(x %in% tetraploid.tumors.discovery){
    "4"
  }else{
    "Other"
  }
})), 
Telomere = telomere.classification.discovery[tumors.discovery],
Sample = sample.information.discovery[tumors.discovery, ]$Sample.type,
Stage = replace(sample.information.discovery$Stage, sample.information.discovery$Stage==5, "4S"),
RiskGroup = sample.information.discovery[tumors.discovery,]$ManualScore)

column_annotation$Sample[column_annotation$Sample %in% "Metastasis"] <- "Primary"
column_annotation$Sample[column_annotation$Sample %in% "Relapse metastasis"] <- "Relapse"

col = c("black", "black", "firebrick", c(brewer.pal(n=8, name="Paired")))
names(col) <- c("Gain", "Loss", "nonsynonymous SNV", "splicing", "stopgain", "frameshift deletion", "nonframeshift deletion",
                 "DUP", "AMP", "DEL", "SV")

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = "white", col = NA))
  },
  "Gain" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["Gain"], col = NA))
  },
  "Loss" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["Loss"], col = NA))
  },
  "AMP" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  "nonsynonymous SNV" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.4, 
              gp = gpar(fill = col["nonsynonymous SNV"], col = NA))
  },
  "DEL" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["DEL"], col = NA))
  },
  "SV" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.4, 
              gp = gpar(fill = col["SV"], col = NA))
  },
  "nonframeshift deletion" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["nonframeshift deletion"], col = NA))
  },
  "splicing" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["splicing"], col = NA))
  },
  "stopgain" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["stopgain"], col = NA))
  },
  "frameshift deletion" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["frameshift deletion"], col = NA))
  },
  "DUP" = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.8, 
              gp = gpar(fill = col["DUP"], col = NA))
  }
)

column_title = "OncoPrint for Neuroblastoma"

telomere.colors <- c(telomere.colors, IR=telomere.colors["IMR"])

annotation_colors <- list(Ploidy = ploidy.cols, Telomere = telomere.colors, 
                          MYCN =  c(`1`="black", `0`="white"), ALT =  c(`1`="black", `0`="white"), TERT =  c(`1`="black", `0`="white"),
                          Sample = time.colors, Stage=stage.colors, Location = time.colors)

mat[mat=="0"] <- ""

mat <- mat[order(rowSums(mat!=""),decreasing=T),]
## split the columns by: HR/LR, Primary, Metastasis, Relapse, Timing of ECA possible and Ploidy

split.criterion <- column_annotation$RiskGroup
split.criterion <- factor(split.criterion, levels=c("LR", "IR", "HR"))

## one oncoprint for tumors with ECA only

row.split.criterion <- ifelse(rownames(mat) %in% c("7 gain", "7q gain", "1q gain", "1 gain", "2 gain", "2p gain", "17 gain", "17q gain", 
                                                   "11q loss", "1p loss"), "CNVs", "SNVs")

mat <- rbind(mat[c("17q gain", "17 gain", "1p loss", "1q gain", "1 gain", "7q gain", "7 gain", "2p gain", "2 gain", "11q loss"), ],
             mat[row.split.criterion=="SNVs",])

addWorksheet(wb, "d")
writeData(wb, "d", mat, rowNames=T)

row.split.criterion <- c(row.split.criterion[row.split.criterion=="CNVs"], row.split.criterion[row.split.criterion=="SNVs"])

oncoPrint(mat, column_split = split.criterion,
          row_split = factor(row.split.criterion, levels=c("CNVs", "SNVs")), row_order = 1:nrow(mat),
          bottom_annotation = HeatmapAnnotation(df=column_annotation, col=annotation_colors),
          alter_fun = alter_fun, col = col, row_names_gp = gpar(fontsize = 8),
          remove_empty_columns = F, remove_empty_rows = TRUE, show_column_names = T,
          column_title = column_title, show_pct = F)

dev.off()



##############################################################################################################################################
## Fig. 1e: Mutational signatures

pdf(paste0(panel.directory, "Figure_1_e.pdf"), width=5, height = 3, useDingbats = F)

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

to.plot <- rbind(sigs.discovery, clonal.sigs.discovery, subclonal.sigs.discovery)
to.plot$Clonality <- c(rep("All", nrow(clonal.sigs.discovery)), rep("Clonal", nrow(clonal.sigs.discovery)), rep("Subclonal", nrow(subclonal.sigs.discovery)))
to.plot <- melt(to.plot, id.vars = c( "mutations", "Sample", "Clonality"), variable.name = "Signature", value.name = "Exposure")
to.plot <- to.plot[,-which(colnames(to.plot)=="mutations")]
to.plot <- reshape(to.plot, timevar = "Clonality", idvar = c( "Sample", "Signature"), direction = "wide")

to.plot$Sample <- factor(to.plot$Sample, levels=colnames(mat))

## stratify by the groups shown in the oncoprint, do for all SNVs only  
p <- list()

p[[1]] <- ggplot(to.plot[!to.plot$Signature %in% c("SBS1", "SBS5", "SBS40") &
                       to.plot$Sample %in% colnames(mat)[split.criterion=="LR"],], aes(x=Sample, y=Exposure.All, fill=Signature)) + 
  geom_col() + scale_fill_manual(values = replace(sig.colors, names(sig.colors)=="Clock", "darkgrey"))

p[[2]] <-  ggplot(to.plot[!to.plot$Signature %in% c("SBS1", "SBS5", "SBS40") &
                        to.plot$Sample %in% colnames(mat)[split.criterion=="IR"],], aes(x=Sample, y=Exposure.All, fill=Signature)) + 
  geom_col() + scale_fill_manual(values = replace(sig.colors, names(sig.colors)=="Clock", "darkgrey")) 

p[[3]] <-  ggplot(to.plot[!to.plot$Signature %in% c("SBS1", "SBS5", "SBS40") &
                        to.plot$Sample %in% colnames(mat)[split.criterion=="HR"],], aes(x=Sample, y=Exposure.All, fill=Signature)) + 
  geom_col() + scale_fill_manual(values = replace(sig.colors, names(sig.colors)=="Clock", "darkgrey")) 


## remove axes
p <- lapply(p, function(x){x +  theme(axis.line=element_blank(),
                                      axis.text.x=element_blank(),
                                      axis.text.y=element_blank(),
                                      axis.ticks=element_blank(),
                                      axis.title.x=element_blank(),
                                      axis.title.y=element_blank())})


ggarrange(plotlist = p, nrow=3, ncol=3, widths = rep(table(split.criterion), 3), common.legend = T)

addWorksheet(wb, "e")
writeData(wb, "e", to.plot) 


dev.off()

##############################################################################################################################################
## Fig. S1e: Mutational signatures correlation between clonal, subclonal, all

pdf(paste0(panel.directory, "Figure_S1e.pdf"), width=6, height = 3, useDingbats = F)

## do for each sig
p <- list()
for(i in unique(to.plot$Signature)){
  if(i %in% c("SBS1", "SBS5", "SBS40")){next} ## summarize these into 1 clock-like sig.
  p[[length(p) + 1]] <- ggplot(to.plot[to.plot$Signature==i,], aes(x=Exposure.All*100, y=Exposure.Clonal*100)) + 
    geom_point(col=sig.colors[i], shape=1, size=0.5) +
    geom_abline(slope = 1, intercept = 0) + ggtitle(i) + scale_x_continuous("% all SNVs", breaks = c(0, 50, 100), limits = c(0,100)) + 
    scale_y_continuous("% clonal SNVs", breaks = c(0, 50, 100), limits = c(0,100)) + stat_cor(size=1, p.accuracy = 0.001) + 
    theme(aspect.ratio = 1)
  
  p[[length(p) + 1]] <- ggplot(to.plot[to.plot$Signature==i,], aes(x=Exposure.All*100, y=Exposure.Subclonal*100)) + 
    geom_point(col=sig.colors[i], shape=1, size=0.5) +
    geom_abline(slope = 1, intercept = 0) + ggtitle(i)+ scale_x_continuous("% all SNVs", breaks = c(0, 50, 100), limits = c(0,100)) + 
    scale_y_continuous("% subclonal SNVs", breaks = c(0, 50, 100), limits = c(0,100))   + stat_cor(size=1, p.accuracy = 0.001) + 
    theme(aspect.ratio = 1)
  
  p[[length(p) + 1]] <- ggplot(to.plot[to.plot$Signature==i,], aes(x=Exposure.Clonal*100, y=Exposure.Subclonal*100)) + 
    geom_point(col=sig.colors[i], shape=1, size=0.5) +
    geom_abline(slope = 1, intercept = 0) + ggtitle(i)+ scale_x_continuous("% clonal SNVs", breaks = c(0, 50, 100), limits = c(0,100)) + 
    scale_y_continuous("% subclonal SNVs", breaks = c(0, 50, 100), limits = c(0,100)) + stat_cor(size=1, p.accuracy = 0.001) + 
    theme(aspect.ratio = 1)
}


ggarrange(plotlist=p, nrow=3, ncol=6)


addWorksheet(wb.s, "e")
writeData(wb.s, "e", to.plot) 


dev.off()


##############################################################################################################################################
saveWorkbook(wb, file = paste0(panel.directory, "Source_data_Fig.1.xlsx"), overwrite=T)
saveWorkbook(wb.s, file = paste0(panel.directory, "Source_data_Fig.S1.xlsx"), overwrite=T)
