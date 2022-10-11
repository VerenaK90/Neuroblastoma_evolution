## Reproduce Fig. 5
##############################################################################################################################################
## load settings and libraries
source("./Settings.R")

load(paste0(rdata.directory, "MRCA_timing.RData"))

## source data:
wb <- createWorkbook()

## store figure panels
panel.directory <- paste0(output.directory, "Figure5/")

if(!dir.exists(panel.directory)){
  dir.create(panel.directory)
}

source(paste0(custom.script.directory, "Survival_analysis.R"))


##########################################################################################################################################
## Figure 5a: Odds ratio for covariates with significant enrichment in early or late MRCAs

to.plot <- joined.p.value.table
to.plot <- to.plot[to.plot$p.value < 0.05,]
to.plot$Parameter <- factor(to.plot$Parameter, levels=c("Sex", "TMM.binary", "Stage.binary", "Age.binary", "Triploidy",
                                                        "7.gain", "17.gain",
                                                        "1p.deletion", "1q.gain", "7q.gain", "11q.deletion", "17q.gain"))

addWorksheet(wb, "5a")
writeData(wb, "5a", to.plot)


pdf(paste0(panel.directory,"Figure_5a.pdf"), useDingbats = F)

ggplot(to.plot, aes(y=odds.ratio, ymin=odds.ratio.l, ymax=odds.ratio.u, x=Parameter)) + coord_flip() + 
  geom_pointrange() + geom_hline(yintercept = 1, linetype=2) + scale_y_log10()

dev.off()

##########################################################################################################################################
## Figure 5b, enrichment of TMM mechanisms among tumors with early and late MRCA

pdf(paste0(panel.directory,"Figure_5b.pdf"), useDingbats = F)

to.plot <- joined.categorized.by.MRCA[joined.categorized.by.MRCA$MRCA.time=="low",]

to.plot <- to.plot %>% group_by(`Telomere.maintenance.mechanism`) %>% count(`Telomere.maintenance.mechanism`)

addWorksheet(wb, "b_early_MRCA")
writeData(wb, "b_early_MRCA", to.plot)


ggplot(to.plot, aes(x="", y=n, fill=Telomere.maintenance.mechanism)) +
  geom_col( width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values=telomere.colors) + ggtitle("Early MRCA")

to.plot <- joined.categorized.by.MRCA[joined.categorized.by.MRCA$MRCA.time=="high",]

to.plot <- to.plot %>% group_by(`Telomere.maintenance.mechanism`) %>% count(`Telomere.maintenance.mechanism`)

addWorksheet(wb, "b_late_MRCA")
writeData(wb, "b_late_MRCA", to.plot)

ggplot(to.plot, aes(x="", y=n, fill=Telomere.maintenance.mechanism)) +
  geom_col( width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values=telomere.colors) + ggtitle("Late MRCA")

dev.off()


##########################################################################################################################################
## Figure 5c, stratify MRCA by TMM

pdf(paste0(panel.directory,"Figure_5c.pdf"), useDingbats = F)

to.plot <- data.frame(MRCA=mutation.time.mrca[primary.tumors.discovery,]$Mean,
                      MRCA.lower = mutation.time.mrca[primary.tumors.discovery,]$Min,
                      MRCA.upper = mutation.time.mrca[primary.tumors.discovery,]$Max,
                      TMM = sample.information.discovery[primary.tumors.discovery,]$Telomere.maintenance.mechanism)

primary.tumors.validation <- rownames(sample.information.validation[sample.information.validation$Sample.type %in% c("Primary", "Metastasis"),])
to.plot <- rbind(to.plot,
                 data.frame(MRCA=mutation.time.mrca[primary.tumors.validation,]$Mean,
                            MRCA.lower = mutation.time.mrca[primary.tumors.validation,]$Min,
                            MRCA.upper = mutation.time.mrca[primary.tumors.validation,]$Max,
                            TMM = sample.information.validation[sample.information.validation$Sample.type %in% c("Primary", "Metastasis"),]$Telomere.maintenance.mechanism) )

to.plot <- to.plot[order(to.plot$MRCA),]

to.plot$MRCA <- to.plot$MRCA/3.3/10^3
to.plot$MRCA.lower <- to.plot$MRCA.lower/3.3/10^3
to.plot$MRCA.upper <- to.plot$MRCA.upper/3.3/10^3

addWorksheet(wb, "c")
writeData(wb, "c", to.plot)


p <- ggplot(data = to.plot[to.plot$TMM=="MNA",],
            aes(x=MRCA, y=seq(1/length(MRCA),1,length.out = length(MRCA)),
                ymin =  sapply(sort(MRCA), function(x){
                  sum(MRCA.upper <= x)
                })/length(MRCA),
                ymax= sapply(sort(MRCA), function(x){
                  sum(MRCA.lower <= x)
                })/length(MRCA),
                fill=TMM
            )) +
  geom_stepribbon(col=NA)+
  stat_ecdf(col="darkgrey") +
  scale_color_manual(values=telomere.colors) + 
  scale_fill_manual(values=alpha(telomere.colors)) +
  scale_x_continuous(name = "SSNVs/Mb")+
  theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(name = "Fraction of tumors") 


p <- p +
  geom_stepribbon(data = to.plot[to.plot$TMM=="ALT",],
                  aes(x=MRCA, y=seq(1/length(MRCA),1,length.out = length(MRCA)),
                      ymin =  sapply(sort(MRCA), function(x){
                        sum(MRCA.upper <= x)
                      })/length(MRCA),
                      ymax= sapply(sort(MRCA), function(x){
                        sum(MRCA.lower <= x)
                      })/length(MRCA)),col=NA)+ 
  stat_ecdf(data = to.plot[to.plot$TMM=="ALT",], aes(x=MRCA, y=seq(1/length(MRCA),1,length.out = length(MRCA))),
            col="darkgrey") 


p <- p + geom_stepribbon(data = to.plot[to.plot$TMM=="TERT",],
                         aes(x=MRCA, y=seq(1/length(MRCA),1,length.out = length(MRCA)),
                             ymin =  sapply(sort(MRCA), function(x){
                               sum(MRCA.upper <= x)
                             })/length(MRCA),
                             ymax= sapply(sort(MRCA), function(x){
                               sum(MRCA.lower <= x)
                             })/length(MRCA)),col="lightgrey") + 
  stat_ecdf(data = to.plot[to.plot$TMM=="TERT",],
            aes(x=MRCA, y=seq(1/length(MRCA),1,length.out = length(MRCA))
            ), col="darkgrey") 



p <- p+ geom_stepribbon(data = to.plot[to.plot$TMM=="None",],
                        aes(x=MRCA, y=seq(1/length(MRCA),1,length.out = length(MRCA)),
                            ymin =  sapply(sort(MRCA), function(x){
                              sum(MRCA.upper <= x)
                            })/length(MRCA),
                            ymax= sapply(sort(MRCA), function(x){
                              sum(MRCA.lower <= x)
                            })/length(MRCA)),col="lightgrey") + stat_ecdf(data = to.plot[to.plot$TMM=="None",],
                                                                          aes(x=MRCA, y=seq(1/length(MRCA),1,length.out = length(MRCA))
                                                                          ),
                                                                          col="darkgrey")



print(p)


dev.off()


##########################################################################################################################################

saveWorkbook(wb, file = paste0(panel.directory,"Source_data_Fig.5.xlsx"), overwrite=T)
