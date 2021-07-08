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
## Figure 4a, stratify MRCA by TMM
source(paste0(custom.script.directory, "Survival_analysis.R"))

pdf(paste0(panel.directory,"Figure_4a.pdf"), useDingbats = F)

to.plot <- data.frame(MRCA=mutation.time.mrca[primary.tumors.80x],
                      MRCA.lower = mutation.time.mrca.lower[primary.tumors.80x],
                      MRCA.upper = mutation.time.mrca.upper[primary.tumors.80x],
                      TMM = sample.information.80x[primary.tumors.80x,]$Telomere.maintenance.mechanism)

primary.tumors.30x <- rownames(sample.information.30x[sample.information.30x$Location %in% c("Primary", "Metastasis"),])
to.plot <- rbind(to.plot,
                 data.frame(MRCA=mutation.time.mrca[primary.tumors.30x],
                            MRCA.lower = mutation.time.mrca.lower[primary.tumors.30x],
                            MRCA.upper = mutation.time.mrca.upper[primary.tumors.30x],
                            TMM = sample.information.30x[sample.information.30x$Location %in% c("Primary", "Metastasis"),]$Telomere.maintenance.mechanism) )

to.plot <- to.plot[order(to.plot$MRCA),]

to.plot$MRCA <- to.plot$MRCA/3.3/10^3
to.plot$MRCA.lower <- to.plot$MRCA.lower/3.3/10^3
to.plot$MRCA.upper <- to.plot$MRCA.upper/3.3/10^3

addWorksheet(wb, "a")
writeData(wb, "a", to.plot)


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
## Figure 4b, enrichment of TMM mechanisms among tumors with early and late MRCA

pdf(paste0(panel.directory,"Figure_4b.pdf"), useDingbats = F)

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

saveWorkbook(wb, file = paste0(panel.directory,"Source_data_Fig.4.xlsx"), overwrite=T)
