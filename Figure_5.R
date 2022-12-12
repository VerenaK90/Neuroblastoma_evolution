## Reproduce Fig. 5
##############################################################################################################################################
## load settings and libraries
source("./Nextcloud/NB_manuscript/Submission_NG/RevisionII//Plots_and_scripts/Custom_scripts/Settings.R")

load(paste0(rdata.directory, "MRCA_timing.RData"))

## source data:
wb <- createWorkbook()
wb.s <- createWorkbook()

## store figure panels
panel.directory <- paste0(output.directory, "Figure5/")

if(!dir.exists(panel.directory)){
  dir.create(panel.directory)
}


source(paste0(custom.script.directory, "Survival_analysis.R"))
##########################################################################################################################################
## Figure 5 a,b: Survival curves, overall survival for validation cohort only

addWorksheet(wb, "a")
writeData(wb, "a",  data.frame(Time=EFS.fit.validation$time, n.Risk = EFS.fit.validation$n.risk, 
                               Category=c(rep(names(EFS.fit.validation$strata)[1], EFS.fit.validation$strata[1]),
                                          rep(names(EFS.fit.validation$strata)[2], EFS.fit.validation$strata[2])),
                               Censored=EFS.fit.validation$n.censor, Event = EFS.fit.validation$n.event, Survival=EFS.fit.validation$surv,
                               LowerCI=EFS.fit.validation$lower, UpperCI=EFS.fit.validation$upper))

chars <- capture.output(EFS.fit.validation.stats)

writeData(wb, sheet = "a", chars, startCol = 10)

addWorksheet(wb, "b")
writeData(wb, "b", data.frame(Time=survival.fit.validation$time, n.Risk = survival.fit.validation$n.risk, 
                              Category=c(rep(names(survival.fit.validation$strata)[1], survival.fit.validation$strata[1]),
                                         rep(names(survival.fit.validation$strata)[2], survival.fit.validation$strata[2])),
                              Censored=survival.fit.validation$n.censor, Event = survival.fit.validation$n.event, Survival=survival.fit.validation$surv,
                              LowerCI=survival.fit.validation$lower, UpperCI=survival.fit.validation$upper))

chars <- capture.output(survival.fit.validation.stats)

writeData(wb, sheet = "b", chars, startCol = 10)


pdf(paste0(panel.directory,"Figure_5a_b.pdf"), useDingbats = F)

ggsurvplot(EFS.fit.validation, data = categorized.by.MRCA.validation, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0,10)) 

ggsurvplot(survival.fit.validation, data = categorized.by.MRCA.validation, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0,10)) 

dev.off()


##########################################################################################################################################
## Figure 5c, d: Survival curves, overall survival, event-free survival

addWorksheet(wb, "c")
writeData(wb, "c", data.frame(Time=joined.EFS.fit$time, n.Risk = joined.EFS.fit$n.risk, 
                              Category=c(rep(names(joined.EFS.fit$strata)[1], joined.EFS.fit$strata[1]),
                                         rep(names(joined.EFS.fit$strata)[2], joined.EFS.fit$strata[2])),
                              Censored=joined.EFS.fit$n.censor, Event = joined.EFS.fit$n.event, Survival=joined.EFS.fit$surv,
                              LowerCI=joined.EFS.fit$lower, UpperCI=joined.EFS.fit$upper))

chars <- capture.output(joined.EFS.fit.stats)

writeData(wb, sheet = "c", chars, startCol = 10)

pdf(paste0(panel.directory,"Figure_5c.pdf"), useDingbats = F)

ggsurvplot(joined.EFS.fit, data = joined.categorized.by.MRCA, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0, 10), break.x.by=5, cex=8, linewidth=0.8) 

dev.off()


addWorksheet(wb, "d")
writeData(wb, "d", data.frame(Time=joined.survival.fit$time, n.Risk = joined.survival.fit$n.risk, 
                              Category=c(rep(names(joined.survival.fit$strata)[1], joined.survival.fit$strata[1]),
                                         rep(names(joined.survival.fit$strata)[2], joined.survival.fit$strata[2])),
                              Censored=joined.survival.fit$n.censor, Event = joined.survival.fit$n.event, Survival=joined.survival.fit$surv,
                              LowerCI=joined.survival.fit$lower, UpperCI=joined.survival.fit$upper))

chars <- capture.output(joined.survival.fit.stats)

writeData(wb, sheet = "d", chars, startCol = 10)

pdf(paste0(panel.directory,"Figure_5d.pdf"), useDingbats = F)
ggsurvplot(joined.survival.fit, data = joined.categorized.by.MRCA, risk.table = TRUE, pval=T, conf.int = T, color="strata", censor.shape=124,
           palette=c("dodgerblue", "dodgerblue4"), xlim=c(0, 10), break.x.by=5, cex=8, linewidth=0.8) 
dev.off()



##############################################################################################################################################
## Figure 4e, S4a: Cox regression

pdf(paste0(panel.directory,"Figure_5e.pdf"), useDingbats = F, width = 4, height=4)

ggforest(fit.coxph_MRCA_TMM_Stage_Age_RAS.EFS, data = joined.categorized.by.MRCA, main = "EFS") 

dev.off()

chars <- capture.output(summary(fit.coxph_MRCA_TMM_Stage_Age_RAS.EFS))

addWorksheet(wb, "e")
writeData(wb, sheet = "e", chars)


pdf(paste0(panel.directory,"Figure_S5a.pdf"), useDingbats = F, width = 4, height=4)

ggforest(fit.coxph_MRCA_TMM_Stage_Age_RAS.OS, data = joined.categorized.by.MRCA, main="OS")

dev.off()


chars <- capture.output(summary(fit.coxph_MRCA_TMM_Stage_Age_RAS.OS))

addWorksheet(wb.s, "a")
writeData(wb.s, sheet = "a", chars)

##########################################################################################################################################

saveWorkbook(wb, file = paste0(panel.directory,"Source_data_Fig.5.xlsx"), overwrite=T)
saveWorkbook(wb.s, file = paste0(panel.directory,"Source_data_Fig.S4.xlsx"), overwrite=T)
