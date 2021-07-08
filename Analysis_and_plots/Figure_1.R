## Reproduce Fig. 1 - detection data set only
##############################################################################################################################################
## load settings and libraries

source(paste0(custom.script.directory, "Gains_and_losses.R"))
source(paste0(custom.script.directory, "Driver_mutations.R"))
source(paste0(custom.script.directory, "Oncoprint_Fig1.R"))
source(paste0(custom.script.directory, "Mutational_signatures.R"))

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

to.plot <- sample.information.80x[,c("Telomere.maintenance.mechanism", "Location", "Age", "ManualScore", "Stage")]
to.plot$Patient_ID <- rownames(to.plot)
## manipulate data slightly:
to.plot$Age.binary <- "<18 months"
to.plot$Age.binary[to.plot$Age >= 18*30]  <- ">=18 months"
to.plot$Age <- to.plot$Age.binary
to.plot$Age <- factor(to.plot$Age, levels=c("<18 months", ">=18 months"))

to.plot$Stage[to.plot$Stage==5] <- "4S"
to.plot$Stage <- factor(to.plot$Stage, levels=c("1", "4S", "2", "3", "4"))

##

to.plot$Location[to.plot$Location=="Relapse 3"] <- "Relapse"
to.plot$Location[to.plot$Location %in% c("Relapse tumor", "Relapse metastasis")] <- "Relapse"
to.plot$Location[to.plot$Location %in% c("Primary", "Metastasis")] <- "Primary"
to.plot$Location <- factor(to.plot$Location, levels=c("Primary", "Relapse"))
to.plot$Telomere.maintenance.mechanism <- factor(to.plot$Telomere.maintenance.mechanism, levels=c("MNA", "TERT", "ALT", "Multiple", "None"))
to.plot$Sample_pair <- "No pair"
to.plot$Sample_pair[to.plot$Patient_ID %in% c("NBE51", "NBE78")] <- "NBE51/NBE78"
to.plot$Sample_pair[to.plot$Patient_ID %in% c("NBE51", "NBE66")] <- "NBE11/NBE66"
to.plot$Sample_pair <- factor(to.plot$Sample_pair)


to.plot <- to.plot[order(to.plot$Telomere.maintenance.mechanism, decreasing = F),]
to.plot <- to.plot[order(to.plot$Age, decreasing = F),]
to.plot <- to.plot[order(to.plot$Location, decreasing = F),]
to.plot <- to.plot[order(to.plot$Stage, decreasing = F),]

addWorksheet(wb, "a")
writeData(wb, "a", to.plot)

pdf(paste0(panel.directory, "/Figure_1a.pdf"))

##separate heatmaps for high-risk, intermediate-risk, low-risk and the 2 sample pairs

ht_list <-  Heatmap(to.plot$Stage, name="Stage", col = stage.colors, width=unit(0.5, "cm"), 
                    split=to.plot[,c("Sample_pair", "ManualScore")], rect_gp = gpar(col="white", lwd=2)) +
  Heatmap(to.plot$Location, name="Sample", col = time.colors, width=unit(0.5, "cm"),  rect_gp = gpar(col="white", lwd=2)) +
  Heatmap(as.numeric(to.plot$Age), name="Age", col = c("lightgrey", "black"),
          heatmap_legend_param = list(labels=levels(to.plot$Age), at=c(1,2)), width=unit(0.5, "cm"), rect_gp = gpar(col="white", lwd=2)) +
  Heatmap(to.plot$Telomere.maintenance.mechanism, name="TMM", col = telomere.colors, width=unit(0.5, "cm"),  rect_gp = gpar(col="white", lwd=2)) 

draw(ht_list, ht_gap=unit(0, "cm"))


dev.off()

##############################################################################################################################################
## Figure 1b: Purity/Ploidy overview

pdf(paste0(panel.directory, "/Figure_1b.pdf"), useDingbats = F, width=4, height=4)

addWorksheet(wb, "b")
writeData(wb, "b", sample.information.80x[,c("Ploidy", "Purity")])

ggplot(sample.information.80x, aes(x=Ploidy, y=Purity, group=Ploidy)) + geom_beeswarm() + scale_y_continuous(limits=c(0,1)) + 
  geom_boxplot(fill=NA)

dev.off()


##############################################################################################################################################
## Figure S1a: Ploidy vs. DNA-index*2

pdf(paste0(panel.directory, "/Figure_S1a.pdf"), useDingbats = F, width=4, height=4)

load(paste0(rdata.directory, "Clonal_mutations_different_ploidies.RData"))

to.plot <- data.frame(DNAIndex = as.numeric(sample.information.80x[tumors.80x,]$`DNA-INDEX`),
                      ACEseq = ploidy[tumors.80x])

addWorksheet(wb.s, "a")
writeData(wb.s, "a", to.plot)

ggplot(to.plot, aes(x=ACEseq, y=DNAIndex*2)) + geom_point() + scale_y_continuous(limits=c(1,max(to.plot)*1.1)) +
  scale_x_continuous(limits=c(1, max(to.plot)*1.1)) + geom_abline(slope=1, intercept = 0, linetype=2)

cor.test(to.plot$DNAIndex*2, to.plot$ACEseq )

dev.off()


##############################################################################################################################################
## Figure S1b: Plot association between ploidy and stage

pdf(paste0(panel.directory, "Figure_S1b.pdf"))

to.plot <- data.frame(Stage = sample.information.80x[tumors.80x,]$Stage,
                      Ploidy = sample.information.80x$Ploidy )
to.plot$Ploidy <- factor(to.plot$Ploidy, levels = c("2", "3", "4"))
to.plot$Stage <-  replace(to.plot$Stage, to.plot$Stage=="5", "4S")
to.plot$Stage <- factor(to.plot$Stage, levels = c("1", "4S", "2", "3", "4"))

addWorksheet(wb.s, "b")
writeData(wb.s, "b", to.plot)


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
## Figure S1c: Gains and losses across the cohort

pdf(paste0(panel.directory, "Figure_S1c.pdf"))

to.store <- rbind(data.frame(Chrom=seqnames(gain.granges), Start = start(gain.granges), End = end(gain.granges),
                       Number.of.cases = gain.granges$Coverage, Type = gain.granges$Type),
                  data.frame(Chrom=seqnames(loss.granges), Start = start(loss.granges), End = end(loss.granges),
                             Number.of.cases = loss.granges$Coverage, Type = loss.granges$Type))

addWorksheet(wb.s, "c")
writeData(wb.s, "c", to.plot) 

plotGrandLinear(c(gain.granges, loss.granges), coord="genome", geom="bar", 
                ymax=c(gain.granges, loss.granges)$Coverage/length(tumors.80x)*100, aes(y=Coverage/length(tumors.80x)*100,
                                                                                    col=Type, fill=Type)) +  theme(panel.background = element_rect(fill = "white", color="black"), panel.grid.major = element_blank(),
                                                                                                                   panel.grid.minor = element_blank()) + geom_hline(yintercept = 0) +
  scale_color_manual(values = c("loss" = "#DA1F28", "gain"= "#B0D8A2")) +scale_fill_manual(values = c("loss" = "firebrick", "gain"= "darkgreen"))+
  scale_y_continuous("% Tumors")
dev.off()


##############################################################################################################################################
## Fig. 1d: Oncoprint


pdf(paste0(panel.directory, "/Figure_1d.pdf"), useDingbats = F, width=10, height = 7)

column_annotation <- data.frame(Ploidy = unlist(sapply(tumors.80x, function(x){
  if(x %in% diploid.tumors.80x){
    "2"
  }else if(x %in% triploid.tumors.80x){
    "3"
  }else if(x %in% tetraploid.tumors.80x){
    "4"
  }else{
    "Other"
  }
})), 
Telomere = telomere.classification.80x[tumors.80x],
Sample = sample.information.80x[tumors.80x, ]$Location,
Stage = replace(sample.information.80x$Stage, sample.information.80x$Stage==5, "4S"),
RiskGroup = sample.information.80x[tumors.80x,]$ManualScore)


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
                          Sample = time.colors, Stage=stage.colors)

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
## Fig. 1e-g: Mutational signatures

pdf(paste0(panel.directory, "Figure_1_e_f_g.pdf"), width=5, height = 3, useDingbats = F)

## all mutations
to.plot <- data.frame(Signatures = tmp$ALL[apply(tmp, 1, max)>0.05],
                      Upper = tmp$ALL.min[apply(tmp, 1, max)>0.05],
                      Lower = tmp$ALL.max[apply(tmp, 1, max)>0.05])

ggplot(data = to.plot,
       aes(x = factor(rownames(signature.contribution.c)[apply(tmp, 1, max)>0.05], levels = rownames(signature.contribution.c)[apply(tmp, 1, max)>0.05]), y = Signatures, ymin = Lower, ymax = Upper)) + 
  geom_col() + geom_errorbar() +
  scale_x_discrete(name="Mutational Signature") + scale_y_continuous(name="Exposure") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=12),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + ggtitle("All mutations")

addWorksheet(wb, "e")
writeData(wb, "e", to.plot) 

## subclonal mutations
to.plot <- data.frame(Signatures = tmp$SC[apply(tmp, 1, max)>0.05],
                      Upper = tmp$SC.min[apply(tmp, 1, max)>0.05],
                      Lower = tmp$SC.max[apply(tmp, 1, max)>0.05])

ggplot(data = to.plot,
       aes(x = factor(rownames(signature.contribution.c)[apply(tmp, 1, max)>0.05], levels = rownames(signature.contribution.sc)[apply(tmp, 1, max)>0.05]), y = Signatures, ymin = Lower, ymax = Upper)) + 
  geom_col() + geom_errorbar() +
  scale_x_discrete(name="Mutational Signature") + scale_y_continuous(name="Exposure") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=12),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + ggtitle("Subclonals")

addWorksheet(wb, "f")
writeData(wb, "f", to.plot) 

## clonal mutations
to.plot <- data.frame(Signatures = tmp$C[apply(tmp, 1, max)>0.05],
                      Upper = tmp$C.min[apply(tmp, 1, max)>0.05],
                      Lower = tmp$C.max[apply(tmp, 1, max)>0.05])

ggplot(data = to.plot,
       aes(x = factor(rownames(signature.contribution.c)[apply(tmp, 1, max)>0.05], levels = rownames(signature.contribution.c)[apply(tmp, 1, max)>0.05]), y = Signatures, ymin = Lower, ymax = Upper)) + 
  geom_col() + geom_errorbar() +
  scale_x_discrete(name="Mutational Signature") + scale_y_continuous(name="Exposure") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=12),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + ggtitle("Clonals")

addWorksheet(wb, "g")
writeData(wb, "g", to.plot) 

dev.off()


saveWorkbook(wb, file = paste0(panel.directory, "Source_data_Fig.1.xlsx"), overwrite=T)
saveWorkbook(wb.s, file = paste0(panel.directory, "Source_data_Fig.S1.xlsx"), overwrite=T)
