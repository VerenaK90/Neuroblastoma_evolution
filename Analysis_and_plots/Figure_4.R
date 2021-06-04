## Reproduce Fig. 4
##############################################################################################################################################
## load settings and libraries

load(paste0(rdata.directory, "MRCA_timing.RData"))

## source data:
wb <- createWorkbook()
wb.s <- createWorkbook()

## store figure panels
panel.directory <- paste0(output.directory, "Figure4/")

if(!dir.exists(panel.directory)){
  dir.create(panel.directory)
}


##########################################################################################################################################
## Figure 4a,b: Survival curves, overall survival, event-free survival

source(paste0(custom.script.directory, "Survival_analysis.R"))

addWorksheet(wb, "a")
writeData(wb, "a", data.frame(Time=joined.survival.fit$time, n.Risk = joined.survival.fit$n.risk, 
                              Category=c(rep(names(joined.survival.fit$strata)[1], joined.survival.fit$strata[1]),
                                         rep(names(joined.survival.fit$strata)[2], joined.survival.fit$strata[2])),
                              Censored=joined.survival.fit$n.censor, Event = joined.survival.fit$n.event, Survival=joined.survival.fit$surv,
                              LowerCI=joined.survival.fit$lower, UpperCI=joined.survival.fit$upper))

pdf(paste0(panel.directory,"Figure_4a.pdf"), useDingbats = F)
ggsurvplot(joined.survival.fit, data = joined.categorized.by.MRCA, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0, 10)) 
dev.off()



addWorksheet(wb, "b")
writeData(wb, "b", data.frame(Time=joined.EFS.fit$time, n.Risk = joined.EFS.fit$n.risk, 
                              Category=c(rep(names(joined.EFS.fit$strata)[1], joined.EFS.fit$strata[1]),
                                         rep(names(joined.EFS.fit$strata)[2], joined.EFS.fit$strata[2])),
                              Censored=joined.EFS.fit$n.censor, Event = joined.EFS.fit$n.event, Survival=joined.EFS.fit$surv,
                              LowerCI=joined.EFS.fit$lower, UpperCI=joined.EFS.fit$upper))

pdf(paste0(panel.directory,"Figure_4b.pdf"), useDingbats = F)

ggsurvplot(joined.EFS.fit, data = joined.categorized.by.MRCA, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0, 10)) 

dev.off()


##########################################################################################################################################
## Figure S3a: Survival curves, overall survival for discovery cohort only


addWorksheet(wb.s, "a_OS")
writeData(wb.s, "a_OS", data.frame(Time=survival.fit$time, n.Risk = survival.fit$n.risk, 
                                   Category=c(rep(names(survival.fit$strata)[1], survival.fit$strata[1]),
                                              rep(names(survival.fit$strata)[2], survival.fit$strata[2])),
                                   Censored=survival.fit$n.censor, Event = survival.fit$n.event, Survival=survival.fit$surv,
                                   LowerCI=survival.fit$lower, UpperCI=survival.fit$upper))

addWorksheet(wb.s, "a_EFS")
writeData(wb.s, "a_EFS",  data.frame(Time=EFS.fit$time, n.Risk = EFS.fit$n.risk, 
                                     Category=c(rep(names(EFS.fit$strata)[1], EFS.fit$strata[1]),
                                                rep(names(EFS.fit$strata)[2], EFS.fit$strata[2])),
                                     Censored=EFS.fit$n.censor, Event = EFS.fit$n.event, Survival=EFS.fit$surv,
                                     LowerCI=EFS.fit$lower, UpperCI=EFS.fit$upper))



pdf(paste0(panel.directory,"Figure_S3a.pdf"), useDingbats = F)
ggsurvplot(survival.fit, data = categorized.by.MRCA, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0,10)) 

ggsurvplot(EFS.fit, data = categorized.by.MRCA, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0,10)) 

dev.off()


##########################################################################################################################################
## Figure S3b: Survival curves, overall survival for validation cohort only


addWorksheet(wb.s, "b_OS")
writeData(wb.s, "b_OS", data.frame(Time=survival.fit.30x$time, n.Risk = survival.fit.30x$n.risk, 
                                    Category=c(rep(names(survival.fit.30x$strata)[1], survival.fit.30x$strata[1]),
                                               rep(names(survival.fit.30x$strata)[2], survival.fit.30x$strata[2])),
                                    Censored=survival.fit.30x$n.censor, Event = survival.fit.30x$n.event, Survival=survival.fit.30x$surv,
                                    LowerCI=survival.fit.30x$lower, UpperCI=survival.fit.30x$upper))

addWorksheet(wb.s, "b_EFS")
writeData(wb.s, "b_EFS",  data.frame(Time=EFS.fit.30x$time, n.Risk = EFS.fit.30x$n.risk, 
                                      Category=c(rep(names(EFS.fit.30x$strata)[1], EFS.fit.30x$strata[1]),
                                                 rep(names(EFS.fit.30x$strata)[2], EFS.fit.30x$strata[2])),
                                      Censored=EFS.fit.30x$n.censor, Event = EFS.fit.30x$n.event, Survival=EFS.fit.30x$surv,
                                      LowerCI=EFS.fit.30x$lower, UpperCI=EFS.fit.30x$upper))



pdf(paste0(panel.directory,"Figure_S3b.pdf"), useDingbats = F)
ggsurvplot(survival.fit.30x, data = categorized.by.MRCA.30x, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0,10)) 

ggsurvplot(EFS.fit.30x, data = categorized.by.MRCA.30x, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0,10)) 

dev.off()



##########################################################################################################################################
## Figure 4c: Determine cutoff for survival times

source(paste0(custom.script.directory, "Survival_analysis.R"))

## illustrate cut point

pdf(paste0(panel.directory,"Figure_4c.pdf"), useDingbats = F)


addWorksheet(wb, "c")
writeData(wb, "c", data.frame(MRCA=joined.categorized.by.MRCA$MRCA))

ggplot(data=joined.categorized.by.MRCA,
       aes(x=MRCA)) + geom_histogram(binwidth=0.02, fill=manual.colors["Late"]) + 
  geom_vline(xintercept = summary(MRCA.cutpoint)[[1]]/3.3/10^3, col="darkgrey", size=1, linetype=2) +
  scale_x_continuous(name="#SSNVs/Mb") + scale_y_continuous(name="# Tumors")

dev.off()


##########################################################################################################################################
## Figure 4d: Illustrate comparability of ECA and early MRCA 


pdf(paste0(panel.directory,"Figure_4d.pdf"), useDingbats = F)

addWorksheet(wb, "d")
writeData(wb, "d", joined.categorized.by.MRCA[,c("MRCA", "ECA", "MRCA.time")])

joined.categorized.by.MRCA$MRCA.time <- factor(joined.categorized.by.MRCA$MRCA.time, levels=c("low", "high"))

ggplot(joined.categorized.by.MRCA, aes(x=paste("MRCA", MRCA.time), y=MRCA)) + geom_boxplot(col=manual.colors["Late"], fill=manual.colors["Late"], alpha=0.5) + 
  geom_beeswarm(col=manual.colors["Late"]) + 
  geom_boxplot(aes(x="ECA", y=ECA),col=manual.colors["Early"],  fill=manual.colors["Early"], alpha=0.5,inherit.aes = F) +
  geom_beeswarm(aes(x="ECA", y=ECA),col=manual.colors["Early"]) + scale_x_discrete(name="", limits=c("ECA", "MRCA low", "MRCA high")) + 
  scale_y_continuous(name="#SSNVs/Mb") 

dev.off()



##########################################################################################################################################
## Figure S3c: Survival based on telomere maintenance mechanism

addWorksheet(wb.s, "c_OS")
writeData(wb.s, "c_OS", data.frame(Time=joined.survival.fit.tmm$time, n.Risk = joined.survival.fit.tmm$n.risk, 
                                   Category=c(rep(names(joined.survival.fit.tmm$strata)[1], joined.survival.fit.tmm$strata[1]),
                                              rep(names(joined.survival.fit.tmm$strata)[2], joined.survival.fit.tmm$strata[2])),
                                   Censored=joined.survival.fit.tmm$n.censor, Event = joined.survival.fit.tmm$n.event, Survival=joined.survival.fit.tmm$surv,
                                   LowerCI=joined.survival.fit.tmm$lower, UpperCI=joined.survival.fit.tmm$upper))

addWorksheet(wb.s, "c_EFS")
writeData(wb.s, "c_EFS", data.frame(Time=joined.EFS.fit.tmm$time, n.Risk = joined.EFS.fit.tmm$n.risk, 
                                    Category=c(rep(names(joined.EFS.fit.tmm$strata)[1], joined.EFS.fit.tmm$strata[1]),
                                               rep(names(joined.EFS.fit.tmm$strata)[2], joined.EFS.fit.tmm$strata[2])),
                                    Censored=joined.EFS.fit.tmm$n.censor, Event = joined.EFS.fit.tmm$n.event, Survival=joined.EFS.fit.tmm$surv,
                                    LowerCI=joined.EFS.fit.tmm$lower, UpperCI=joined.EFS.fit.tmm$upper))


pdf(paste0(panel.directory,"Figure_S3c.pdf"), useDingbats = F)
ggsurvplot(joined.survival.fit.tmm, data = joined.categorized.by.MRCA, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0, 10)) 



ggsurvplot(joined.EFS.fit.tmm, data = joined.categorized.by.MRCA, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0, 10)) 

dev.off()




##########################################################################################################################################
## Figure 4e: Odds ratio for covariates with significant enrichment in early or late MRCAs

to.plot <- joined.p.value.table
to.plot <- to.plot[to.plot$p.value < 0.05,]
to.plot$Parameter <- factor(to.plot$Parameter, levels=c("Sex", "TMM.binary", "Stage.binary", "Age.binary", "Triploidy",
                                                        "1p.deletion", "1q.gain", "7q.gain", "7.gain", "11q.deletion", "17q.gain", "17.gain"))

addWorksheet(wb, "e")
writeData(wb, "e", to.plot)


pdf(paste0(panel.directory,"Figure_4e.pdf"), useDingbats = F)

ggplot(to.plot, aes(y=odds.ratio, ymin=odds.ratio.l, ymax=odds.ratio.u, x=Parameter)) + coord_flip() + 
  geom_pointrange() + geom_hline(yintercept = 1, linetype=2) + scale_y_log10()

dev.off()


##########################################################################################################################################
## Figure 4f, enrichment of TMM mechanisms among tumors with early and late MRCA

pdf(paste0(panel.directory,"Figure_4f.pdf"), useDingbats = F)

to.plot <- joined.categorized.by.MRCA[joined.categorized.by.MRCA$MRCA.time=="low",]

to.plot <- to.plot %>% group_by(`Telomere.maintenance.mechanism`) %>% count(`Telomere.maintenance.mechanism`)

addWorksheet(wb, "f_early_MRCA")
writeData(wb, "f_early_MRCA", to.plot)


ggplot(to.plot, aes(x="", y=n, fill=Telomere.maintenance.mechanism)) +
  geom_col( width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values=telomere.colors) + ggtitle("Early MRCA")

to.plot <- joined.categorized.by.MRCA[joined.categorized.by.MRCA$MRCA.time=="high",]

to.plot <- to.plot %>% group_by(`Telomere.maintenance.mechanism`) %>% count(`Telomere.maintenance.mechanism`)

addWorksheet(wb, "f_late_MRCA")
writeData(wb, "f_late_MRCA", to.plot)

ggplot(to.plot, aes(x="", y=n, fill=Telomere.maintenance.mechanism)) +
  geom_col( width=1) +
  coord_polar("y", start=0) +
  scale_fill_manual(values=telomere.colors) + ggtitle("Late MRCA")

dev.off()



##########################################################################################################################################
## Figure 4g, stratify MRCA by TMM

pdf(paste0(panel.directory,"Figure_4g.pdf"), useDingbats = F)

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

addWorksheet(wb, "g")
writeData(wb, "g", to.plot)


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

saveWorkbook(wb, file = paste0(panel.directory,"Source_data_Fig.4.xlsx"), overwrite=T)
saveWorkbook(wb.s, file = paste0(panel.directory,"Source_data_Fig.S3.xlsx"), overwrite=T)
