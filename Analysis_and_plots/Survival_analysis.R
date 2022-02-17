##########################################################################################################################################
### Stratify survival by MRCA time

if("ECA" %in% colnames(sample.information.80x)){
  colnames(sample.information.80x)[colnames(sample.information.80x)=="ECA"] <- "ECA.exists" 
}

survival.information.80x <- sample.information.80x[sample.information.80x$Location %in% c("Primary", "Metastasis"),
                          c("OS.time.(days)", "OS", "EFS", "EFS.time(days)", "ECA.exists", "Ploidy", "Telomere.maintenance.mechanism",
                            "17q.gain", "17.gain", "7.gain", "7q.gain", "1.gain", "1q.gain", "1p.deletion", "11q.deletion", "Sex", 
                            "Stage", "Age", "NB2004", "MYCN.(current.disease.episode)")]

load(paste0(rdata.directory, "MRCA_timing.RData"))
mutation.time.eca[names(earliest.mutation.time)] <- earliest.mutation.time
mutation.time.eca.lower[names(earliest.mutation.time)] <- earliest.mutation.time.lower
mutation.time.eca.upper[names(earliest.mutation.time)] <- earliest.mutation.time.upper

rownames(survival.information.80x) <- rownames(sample.information.80x)[sample.information.80x$Location %in% c("Primary", "Metastasis")]
colnames(survival.information.80x) <- c("OS.time.", "OS", "EFS", "EFS.time", "ECA.exists", "Ploidy", "Telomere.maintenance.mechanism",
                                    "17q.gain", "17.gain", "7.gain", "7q.gain", "1.gain", "1q.gain", "1p.deletion", "11q.deletion", "Sex", 
                                    "Stage", "Age", "NB2004", "MYCN.(current.disease.episode)")
survival.information.80x$Ploidy[survival.information.80x$Ploidy %in% c(2,4)] <- "2,4"

survival.information.80x$MRCA.time <- mutation.time.mrca[rownames(survival.information.80x)]
survival.information.80x$ECA.time <- mutation.time.eca[rownames(survival.information.80x)]

### natural cutpoint:
MRCA.cutpoint <- 0.05

##########################################################################################################################################
### Categorize the data accordingly

categorized.by.MRCA <- survival.information.80x
categorized.by.MRCA$`OS.time.` <- categorized.by.MRCA$`OS.time.`/365
categorized.by.MRCA$`EFS.time` <- categorized.by.MRCA$`EFS.time`/365
categorized.by.MRCA$MRCA.time <- ifelse(categorized.by.MRCA$MRCA.time/3.3/10^3 < MRCA.cutpoint, "low", "high")
categorized.by.MRCA$MRCA.time <- factor(categorized.by.MRCA$MRCA.time, levels=c("low", "high"))

survival.fit <- survfit(Surv(`OS.time.`, OS) ~ MRCA.time,
               data = categorized.by.MRCA)

EFS.fit <- survfit(Surv(`EFS.time`, EFS) ~ MRCA.time,
               data = categorized.by.MRCA)


##########################################################################################################################################
### Enrichment test

categorized.by.MRCA$TMM.binary <- ifelse(categorized.by.MRCA$Telomere.maintenance.mechanism=="None", "no TMM", "TMM")
categorized.by.MRCA$Triploidy <- ifelse(categorized.by.MRCA$Ploidy==3, "T", "F")
categorized.by.MRCA$Age.binary <- ifelse(categorized.by.MRCA$Age/30 < 18, "< 18", ">= 18")
categorized.by.MRCA$Stage.binary <- ifelse(categorized.by.MRCA$Stage=="4", "4", "< 4")
p.value.table <- data.frame(Parameter = c("TMM.binary", "17q.gain", "17.gain", "7q.gain", "7.gain", "1p.deletion", "1.gain", "1q.gain", 
                                          "11q.deletion", "Age.binary",  "Stage.binary", "Triploidy"),#, "ALK"), 
                            p.value = Inf,
                            odds.ratio = Inf,
                            odds.ratio.l = Inf,
                            odds.ratio.u = Inf)

for(i in 1:nrow(p.value.table)){
  cont.table <- table(categorized.by.MRCA[,c("MRCA.time", as.character(p.value.table$Parameter[i]))])
  out <- fisher.test(cont.table)
  p.value.table[i,]$p.value <- out$p.value
  p.value.table[i,]$odds.ratio <- out$estimate
  p.value.table[i,]$odds.ratio.l <- out$conf.int[1]
  p.value.table[i,]$odds.ratio.u <- out$conf.int[2]
  
}



##########################################################################################################################################
### 30x dataset

if("ECA" %in% colnames(sample.information.30x)){
  colnames(sample.information.30x)[colnames(sample.information.30x)=="ECA"] <- "ECA.exists" 
}

survival.information.30x <- sample.information.30x[sample.information.30x$Location %in% c("Primary", "Metastasis"),
                                                   c("OS.time.(days)", "OS", "EFS", "EFS.time(days)", "ECA.exists", "Telomere.maintenance.mechanism")]

load(paste0(rdata.directory, "MRCA_timing.RData"))
mutation.time.eca[names(earliest.mutation.time)] <- earliest.mutation.time
mutation.time.eca.lower[names(earliest.mutation.time)] <- earliest.mutation.time.lower
mutation.time.eca.upper[names(earliest.mutation.time)] <- earliest.mutation.time.upper

rownames(survival.information.30x) <- rownames(sample.information.30x)[sample.information.30x$Location %in% c("Primary", "Metastasis")]
colnames(survival.information.30x) <- c("OS.time.", "OS", "EFS", "EFS.time", "ECA.exists", "Telomere.maintenance.mechanism")

survival.information.30x$MRCA.time <- mutation.time.mrca[rownames(survival.information.30x)]
survival.information.30x$ECA.time <- mutation.time.eca[rownames(survival.information.30x)]

##########################################################################################################################################
### Categorize the data 

categorized.by.MRCA.30x <- survival.information.30x
categorized.by.MRCA.30x$`OS.time.` <- categorized.by.MRCA.30x$`OS.time.`/365
categorized.by.MRCA.30x$`EFS.time` <- categorized.by.MRCA.30x$`EFS.time`/365
categorized.by.MRCA.30x$MRCA.time <- ifelse(categorized.by.MRCA.30x$MRCA.time/3.3/10^3 < MRCA.cutpoint, "low", "high")
categorized.by.MRCA.30x$MRCA.time <- factor(categorized.by.MRCA.30x$MRCA.time, levels=c("low", "high"))

survival.fit.30x <- survfit(Surv(`OS.time.`, OS) ~ MRCA.time,
                        data = categorized.by.MRCA.30x)

EFS.fit.30x <- survfit(Surv(`EFS.time`, EFS) ~ MRCA.time,
                   data = categorized.by.MRCA.30x)



##########################################################################################################################################
### Combine datasets

survival.information.30x <- sample.information.30x[sample.information.30x$Location %in% c("Primary", "Metastasis"),]
survival.information.30x$Ploidy <- survival.information.30x$Rounded.ploidy
survival.information.30x$MRCA.time <- mutation.time.mrca[rownames(survival.information.30x)]
survival.information.30x$ECA.time <- mutation.time.eca[rownames(survival.information.30x)]

colnames(survival.information.30x)[colnames(survival.information.30x)=="OS.time.(days)"] <- "OS.time."
colnames(survival.information.30x)[colnames(survival.information.30x)=="EFS.time(days)"] <- "EFS.time"
colnames(survival.information.30x)[colnames(survival.information.30x)=="Age.(days)"] <- "Age"

## determine the optimal cut point of the MRCA 

joined.survival.information <- rbind(survival.information.80x[, c("OS.time.", "OS", "EFS", "EFS.time", "ECA.exists", "Ploidy", "Telomere.maintenance.mechanism",
                                                              "17q.gain", "17.gain", "7.gain", "7q.gain", "1.gain", "1q.gain", "1p.deletion", "11q.deletion", "Sex", 
                                                              "Stage", "Age", "MRCA.time", "NB2004", "MYCN.(current.disease.episode)")], 
                                     survival.information.30x[, c("OS.time.", "OS", "EFS", "EFS.time", "ECA.exists", "Ploidy", "Telomere.maintenance.mechanism",
                                                                         "17q.gain", "17.gain", "7.gain", "7q.gain", "1.gain", "1q.gain", "1p.deletion", "11q.deletion", "Sex", 
                                                                         "Stage", "Age", "MRCA.time", "NB2004", "MYCN.(current.disease.episode)")])



##########################################################################################################################################
### Categorize the data 

joined.categorized.by.MRCA <- joined.survival.information
joined.categorized.by.MRCA$`OS.time.` <- joined.categorized.by.MRCA$`OS.time.`/365
joined.categorized.by.MRCA$`EFS.time` <- joined.categorized.by.MRCA$`EFS.time`/365
joined.categorized.by.MRCA$MRCA.time <- ifelse(joined.categorized.by.MRCA$MRCA.time/3.3/10^3 < MRCA.cutpoint, "low", "high")
joined.categorized.by.MRCA$MRCA.time <- factor(joined.categorized.by.MRCA$MRCA.time, levels=c("low", "high"))


joined.categorized.by.MRCA$TMM.binary <- ifelse(joined.categorized.by.MRCA$Telomere.maintenance.mechanism=="None", "no TMM", "TMM")
joined.categorized.by.MRCA$Triploidy <- ifelse(joined.categorized.by.MRCA$Ploidy==3, "T", "F")
joined.categorized.by.MRCA$Age.binary <- ifelse(joined.categorized.by.MRCA$Age/30 < 18, "< 18", ">= 18")
joined.categorized.by.MRCA$Stage.binary <- ifelse(joined.categorized.by.MRCA$Stage=="4", "4", "< 4")

joined.survival.fit <- survfit(Surv(`OS.time.`, OS) ~ MRCA.time,
                        data = joined.categorized.by.MRCA)

joined.EFS.fit <- survfit(Surv(`EFS.time`, EFS) ~ MRCA.time,
                   data = joined.categorized.by.MRCA)



joined.categorized.by.MRCA$TMM.binary <- ifelse(joined.categorized.by.MRCA$Telomere.maintenance.mechanism == "None", "no TMM", "TMM")
joined.categorized.by.MRCA$MRCA <- c(survival.information.80x$MRCA.time, survival.information.30x$MRCA.time)/3.3/10^3
joined.categorized.by.MRCA$ECA <- c(survival.information.80x$ECA.time, survival.information.30x$ECA.time)/3.3/10^3

joined.p.value.table <- data.frame(Parameter = c("TMM.binary", "17q.gain", "17.gain", "7q.gain", "7.gain", "1p.deletion", "1.gain", "1q.gain", 
                                          "11q.deletion", "Age.binary", "Stage.binary", "Triploidy"),#, "ALK"), 
                            p.value = Inf,
                            odds.ratio = Inf,
                            odds.ratio.l = Inf,
                            odds.ratio.u = Inf)

for(i in 1:nrow(joined.p.value.table)){
  cont.table <- table(joined.categorized.by.MRCA[,c("MRCA.time", as.character(joined.p.value.table$Parameter[i]))])
  out <- fisher.test(cont.table)
  joined.p.value.table[i,]$p.value <- out$p.value
  joined.p.value.table[i,]$odds.ratio <- out$estimate
  joined.p.value.table[i,]$odds.ratio.l <- out$conf.int[1]
  joined.p.value.table[i,]$odds.ratio.u <- out$conf.int[2]
  
}


##########################################################################################################################################

# Fit a Cox proportional hazards model

joined.categorized.by.MRCA$ECA.exists <- as.character(as.logical(joined.categorized.by.MRCA$ECA.exists))
joined.categorized.by.MRCA$NB2004 <- factor(joined.categorized.by.MRCA$NB2004, levels=c("observation", "MRG", "HR"))
joined.categorized.by.MRCA$TMM.binary <- factor(joined.categorized.by.MRCA$TMM.binary, levels=c("no TMM", "TMM"))


## I. univariate Cox regression
covariates <- c("MRCA.time", "TMM.binary", "NB2004.binary")

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(joined.categorized.by.MRCA$`OS.time.`, joined.categorized.by.MRCA$OS)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = joined.categorized.by.MRCA)})

# Extract result 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coefficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)


## II. Multivariate analysis
## check assumptions: 
## - proportional hazards: ratio of hazard should be constant over time between groups

lapply(univ_models, cox.zph)
## non significant for any of the covariates; thus the assumption holds.
for(i in 1:length(univ_models)){
  test.ph <- cox.zph(univ_models[[i]])
  print(ggcoxzph(test.ph))
}

## - linear relationship between covariates and log hazard h(t)

## not an issue for categorial variables and we only have categorial variables

## - influential observations (outliers)

## plots the estimated changes in the regression coefficient after deleting each observation

ggcoxdiagnostics(univ_models[[1]], type = "dfbeta", linear.predictions = F)
ggcoxdiagnostics(univ_models[[2]], type = "dfbeta", linear.predictions = F)
ggcoxdiagnostics(univ_models[[3]], type = "dfbeta", linear.predictions = F)


## Event-free survival
surv_object <- Surv(time=joined.categorized.by.MRCA$`EFS.time`, event=joined.categorized.by.MRCA$EFS) 


fit.coxph.efs <- coxph(surv_object ~ MRCA.time + TMM.binary + NB2004 , 
                       data = joined.categorized.by.MRCA)
