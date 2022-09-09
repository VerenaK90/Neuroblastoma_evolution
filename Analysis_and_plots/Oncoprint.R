##############################################################################################################################################
## read in the driver mutation information
source("Settings.R")

##############################################################################################################################################

## in addition, add chromosomal gains at early and late stages
tumors.with.17q <- rownames(sample.information.discovery)[sample.information.discovery$`17q.gain`=="T"]
tumors.with.whole.17.gain <- rownames(sample.information.discovery)[sample.information.discovery$`17.gain`=="T" & sample.information.discovery$`17q.gain`=="F"]
tumors.with.7q <- rownames(sample.information.discovery)[sample.information.discovery$`7q.gain`=="T"]
tumors.with.whole.7.gain <- rownames(sample.information.discovery)[sample.information.discovery$`7.gain`=="T" & sample.information.discovery$`7q.gain`=="F"]
tumors.with.1p.loss <- rownames(sample.information.discovery)[sample.information.discovery$`1p.deletion`=="T"]
tumors.with.1q.gain <- rownames(sample.information.discovery)[sample.information.discovery$`1q.gain`=="T"]
tumors.with.whole.1.gain <- rownames(sample.information.discovery)[sample.information.discovery$`1.gain`=="T" & sample.information.discovery$`1q.gain`=="F"]
tumors.with.2p <- rownames(sample.information.discovery)[sample.information.discovery$`2p.gain`=="T"]
tumors.with.whole.2.gain <- rownames(sample.information.discovery)[sample.information.discovery$`2.gain`=="T" & sample.information.discovery$`2p.gain`=="F"]
tumors.with.11q.loss <- rownames(sample.information.discovery)[sample.information.discovery$`11q.deletion`=="T"]

not.dateable <- data.frame(Gain = "17q",
                           Sample = rownames(sample.information.discovery)[sample.information.discovery$`17q.gain`=="T" & sample.information.discovery$`17q.gain.<.5`=="F"])
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "17",
                           Sample = rownames(sample.information.discovery)[sample.information.discovery$`17.gain`=="T" & sample.information.discovery$`17.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "7q",
                                 Sample = rownames(sample.information.discovery)[sample.information.discovery$`7q.gain`=="T" & sample.information.discovery$`7q.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "7",
                                 Sample = rownames(sample.information.discovery)[sample.information.discovery$`7.gain`=="T" & sample.information.discovery$`7.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "1q",
                                 Sample = rownames(sample.information.discovery)[sample.information.discovery$`1q.gain`=="T" & sample.information.discovery$`1q.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "1",
                                 Sample = rownames(sample.information.discovery)[sample.information.discovery$`1.gain`=="T" & sample.information.discovery$`1.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "2p",
                                 Sample = rownames(sample.information.discovery)[sample.information.discovery$`2p.gain`=="T" & sample.information.discovery$`2p.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "2",
                                 Sample = rownames(sample.information.discovery)[sample.information.discovery$`2.gain`=="T" & sample.information.discovery$`2.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "11q",
                                 Sample = rownames(sample.information.discovery)[sample.information.discovery$`11q.deletion`=="T"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "1p",
                                 Sample = rownames(sample.information.discovery)[sample.information.discovery$`1p.deletion`=="T"]))


## manual annotation of subclonal gains

subclonal.gains <- data.frame(CHROM=c("17", "7q", "7q", "7q"), SAMPLE=c("NBE18", "NBE77", "NBE17", "NBE47"))


## read in annotation
gains.losses.eca <- data.frame(CHROM=1, SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_whole_1 %in% c("ECA", "ECA/MRCA"),]))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="1p", 
                                                       SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_1p %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="1q", 
                                                         SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_1q %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="2", 
                                                         SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_whole_2 %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="2p", 
                                                         SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_2p %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="7", 
                                                         SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_whole_7 %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="7q", 
                                                         SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_7q %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="11q", 
                                                       SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_11q %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="17", 
                                                         SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_whole_17 %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="17q", 
                                                         SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_17q %in% c("ECA", "ECA/MRCA"),])))

gains.losses.early.tetraploid <- data.frame(CHROM=c(), SAMPLE=c())
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, 
                                       data.frame(CHROM=rep(1, sum(sample.information.discovery$Timing_whole_1 %in% c("early tetraploid"))), SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_whole_1 %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, 
                                       data.frame(CHROM=rep("1p", sum(sample.information.discovery$Timing_1p %in% c("early tetraploid"))), SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_1p %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("1q", sum(sample.information.discovery$Timing_1q %in% c("early tetraploid"))), 
                                                       SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_1q %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("2", sum(sample.information.discovery$Timing_whole_2 %in% c("early tetraploid"))), 
                                                       SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_whole_2 %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("2p", sum(sample.information.discovery$Timing_2p %in% c("early tetraploid"))), 
                                                       SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_2p %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("7", sum(sample.information.discovery$Timing_whole_7 %in% c("early tetraploid"))), 
                                                       SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_whole_7 %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("7q", sum(sample.information.discovery$Timing_7q %in% c("early tetraploid"))), 
                                                       SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_7q %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, 
                                       data.frame(CHROM=rep("11q", sum(sample.information.discovery$Timing_11q %in% c("early tetraploid"))), SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_11q %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("17", sum(sample.information.discovery$Timing_whole_17 %in% c("early tetraploid"))), 
                                                       SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_whole_17 %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("17q", sum(sample.information.discovery$Timing_17q %in% c("early tetraploid"))), 
                                                       SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_17q %in% c("early tetraploid"),])))


gains.losses.mrca <- data.frame(CHROM=1, SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_whole_1 %in% c("MRCA", "ECA/MRCA"),]))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="1p", 
                                                         SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_1p %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="1q", 
                                                         SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_1q %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="2", 
                                                         SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_whole_2 %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="2p", 
                                                         SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_2p %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="7", 
                                                         SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_whole_7 %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="7q", 
                                                         SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_7q %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="11q", 
                                                         SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_11q %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="17", 
                                                         SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_whole_17 %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="17q", 
                                                         SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_17q %in% c("MRCA", "ECA/MRCA"),])))


gains.losses.late.tetraploid <- data.frame(CHROM=c(), SAMPLE=c())
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, 
                                       data.frame(CHROM=rep(1, sum(sample.information.discovery$Timing_whole_1 %in% c("late tetraploid"))), SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_whole_1 %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("1q", sum(sample.information.discovery$Timing_1q %in% c("late tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_1q %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("1p", sum(sample.information.discovery$Timing_1p %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_1p %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("2", sum(sample.information.discovery$Timing_whole_2 %in% c("late tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_whole_2 %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("2p", sum(sample.information.discovery$Timing_2p %in% c("late tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_2p %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("7", sum(sample.information.discovery$Timing_whole_7 %in% c("late tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_whole_7 %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("7q", sum(sample.information.discovery$Timing_7q %in% c("late tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_11q %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("11q", sum(sample.information.discovery$Timing_11q %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_7q %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("17", sum(sample.information.discovery$Timing_whole_17 %in% c("late tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_whole_17 %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("17q", sum(sample.information.discovery$Timing_17q %in% c("late tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.discovery[sample.information.discovery$Timing_17q %in% c("late tetraploid"),])))

##############################################################################################################################################
#### Test whether known driver mutations - in ALK, NRAS, KRAS, HRAS, ROS1 arose before or after the ECA.

prob.early.late.sc <- t(apply(filtered.functional.mutations[filtered.functional.mutations$GENE %in% rownames(mat) &
                                                     !filtered.functional.mutations$EXONIC_CLASSIFICATION %in% c("synonymous SNV", "unknown") &
                                                       filtered.functional.mutations$ANNOVAR_FUNCTION %in% c("exonic", "splicing") ,], 1,
                              function(x){
                                if(is.na(x["TCN"]) | x["TCN"]=="sub"){
                                  return(c(NA, NA, NA, x["SAMPLE"]))
                                }
                                ## probability of that mutation to be early clonal, late clonal, subclonal
                                p.late.clonal <- dbinom(x = as.numeric(x["READS_ALT"]), size =as.numeric(x["READS_REF"]) + as.numeric(x["READS_ALT"]), prob = as.numeric(x["Purity"])/(as.numeric(x["TCN"])*as.numeric(x["Purity"]) + (1-as.numeric(x["Purity"]))*2))
                                
                                if(as.numeric(x["TCN"]) > 2 | (as.numeric(x["TCN"])==2 & as.numeric(x["A"])==2)){
                                  p.early.clonal <- t(sum(dbinom(x = as.numeric(x["READS_ALT"]), size = as.numeric(x["READS_REF"]) + as.numeric(x["READS_ALT"]), prob = (2:as.numeric(x["TCN"]))*as.numeric(x["Purity"])/
                                                                   (as.numeric(x["TCN"])*as.numeric(x["Purity"]) + (1-as.numeric(x["Purity"]))*2))))
                                }else{
                                  p.early.clonal <- p.late.clonal
                                }
                                p.subclonal <- 1-pbinom(q = as.numeric(x["READS_ALT"]), size = as.numeric(x["READS_REF"]) + as.numeric(x["READS_ALT"]), prob = as.numeric(x["Purity"])/
                                                          (as.numeric(x["TCN"])*as.numeric(x["Purity"]) + (1-as.numeric(x["Purity"]))*2))
                                sample <- x["SAMPLE"]
                                
                                return(c(p.early.clonal, p.late.clonal, p.subclonal, sample))
                              }))

colnames(prob.early.late.sc) <- c("Early", "Late", "Subclonal", "Sample")
rownames(prob.early.late.sc) <- (filtered.functional.mutations[filtered.functional.mutations$GENE %in% rownames(mat) &
                                                        !filtered.functional.mutations$EXONIC_CLASSIFICATION %in% c("synonymous SNV", "unknown") &
                                                          filtered.functional.mutations$ANNOVAR_FUNCTION %in% c("exonic", "splicing") ,]$GENE)

prob.early.late.sc <- as.data.frame(prob.early.late.sc)
prob.early.late.sc$Early <- as.numeric(as.character(prob.early.late.sc$Early))
prob.early.late.sc$Late <- as.numeric(as.character(prob.early.late.sc$Late))
prob.early.late.sc$Subclonal <- as.numeric(as.character(prob.early.late.sc$Subclonal))

## select the most likely one for each

prob.early.late.sc$MostLikely <- apply(prob.early.late.sc[,1:3], 1, function(x){
  if(is.na(x[1])){
    "NA"
  }else if(x[3]>0.95){
    "Subclonal"
  }else if(x[2]/(x[1]+x[2]) > 0.5){
    "Early or late clonal" ## could also lie on the non-amplified allele
  }else if(x[2]/(x[1]+x[2]) < 0.5){
    "Early clonal"
  }else{
    "Early or late clonal"
  }
})
## subclonal mutations are these with p(subclonal) > 0.95. 

## Now, for each mutation that is annotated as early clonal, I have to check whether the chromosome it lies on is amplified at the ECA
## or at the MRCA and assign it as either early clonal (if < ECA) or unknown
prob.early.late.sc[prob.early.late.sc$MostLikely=="Early clonal",]
## HRAS on 18533, ALK in 22998 ok, BRAF in 19642 ok, take HRAS of NB-S-570 out, FBN2 also unclear, take out

prob.early.late.sc[c("HRAS.1", "FBN2.2"),"MostLikely"] <- "Early or late clonal"

##############################################################################################################################################
#### design a mutation matrix, where the early/late information is included
matrix.early.late <- matrix("", nrow=10+length(unique(driver.genes)),
                            ncol=length(tumors.discovery), dimnames=list(c(unique(driver.genes),
                                                                 "7", "2", "1", "17", "7q", "2p", "1q", "17q", "11q", "1p"),
                                                               tumors.discovery))

genes <- sapply(rownames(prob.early.late.sc), function(x){strsplit(x, split="[.]")[[1]][1]})

for(i in rownames(matrix.early.late)){
  if(i %in% c("1", "2", "7", "17", "1q", "2p", "7q", "17q", "11q", "1p")){
    matrix.early.late[i,as.character(gains.losses.eca[gains.losses.eca$CHROM==i,]$SAMPLE)] <- "early"
    matrix.early.late[i,as.character(gains.losses.mrca[gains.losses.mrca$CHROM==i,]$SAMPLE)] <- paste(matrix.early.late[i,as.character(gains.losses.mrca[gains.losses.mrca$CHROM==i,]$SAMPLE)], "late", sep=" ")
    matrix.early.late[i,as.character(subclonal.gains[subclonal.gains$CHROM==i,]$SAMPLE)] <- "Subclonal"
    matrix.early.late[i,as.character(gains.losses.early.tetraploid[gains.losses.early.tetraploid$CHROM==i,]$SAMPLE)] <- "early tetraploid"
    matrix.early.late[i,as.character(gains.losses.late.tetraploid[gains.losses.late.tetraploid$CHROM==i,]$SAMPLE)] <- paste(matrix.early.late[i,as.character(gains.losses.late.tetraploid[gains.losses.late.tetraploid$CHROM==i,]$SAMPLE)], "late tetraploid", sep=" ")
    
    
  }else{
    for(j in colnames(matrix.early.late)){
      tmp <- prob.early.late.sc[genes==i & prob.early.late.sc$Sample==j,]
      if(nrow(tmp)==0){
        next
      }
      #tmp <- unlist(tmp)
      matrix.early.late[i,j] <- paste(tmp[,5], collapse = ";")
      
    }
    ## accept amplifications/deletions/translocations as drivers
    matrix.early.late[i,as.character(deletions[deletions$gene==i,]$Sample)] <- paste(matrix.early.late[i,as.character(deletions[deletions$gene==i,]$Sample)], "DEL", "Known_driver", sep=";")
    matrix.early.late[i,as.character(amplifications[amplifications$gene==i,]$Sample)] <- paste(matrix.early.late[i,as.character(amplifications[amplifications$gene==i,]$Sample)], "AMP", "Known_driver", sep=";")
    matrix.early.late[i,as.character(translocations[translocations$gene1==i | translocations$gene2==i,]$Sample)] <- paste(matrix.early.late[i,as.character(translocations[translocations$gene1==i | translocations$gene2==i,]$Sample)], translocations[translocations$gene1==i | translocations$gene2==i,]$SV, "Known_driver", sep=";")
    matrix.early.late[i,as.character(filtered.functional.mutations[filtered.functional.mutations$GENE==i ,]$SAMPLE)] <- paste(matrix.early.late[i,as.character(filtered.functional.mutations[filtered.functional.mutations$GENE==i ,]$SAMPLE)], filtered.functional.mutations[filtered.functional.mutations$GENE==i,]$Known_driver, sep=";")
    
  }
  
}


matrix.early.late[matrix.early.late==0] <- ""

## manual adjustments to include ATRX, MYCN, TERT information from Hartlieb et al., 2020
to.replace <- sample.information.discovery[colnames(matrix.early.late), "ATRX"]## check the nonframeshift deletion
to.replace[to.replace%in% c("SNV", "stopgain", "splicing", "whole chromosome loss", "whole chromosome loss; stopgain",
                            "nonsynonymous SNV")] <- ""
to.replace[to.replace=="DEL; SNV"] <- "DEL"
to.replace[to.replace=="DEL; nonsynonymous SNV"] <- "DEL"
to.replace[to.replace=="DEL; nonframeshift deletion"] <- "DEL"
to.replace[to.replace=="nonframeshift deletion"] <- "DEL"
matrix.early.late["ATRX",] <- paste(matrix.early.late["ATRX",], to.replace, sep=";")
matrix.early.late[is.na(matrix.early.late)] <- ""
matrix.early.late[matrix.early.late %in% c(";wt", ";NA")] <- ""
matrix.early.late[matrix.early.late=="whole chromosome loss; stopgain"] <- "stopgain"
matrix.early.late["MYCN",sample.information.discovery[colnames(matrix.early.late),"Subtype"]=="MNA"] <- paste(matrix.early.late["MYCN",sample.information.discovery[colnames(matrix.early.late),"Subtype"]=="MNA"], "AMP", sep=";")
matrix.early.late["TERT",sample.information.discovery[colnames(matrix.early.late),"Subtype"]=="TERT"] <- paste(matrix.early.late["TERT",sample.information.discovery[colnames(matrix.early.late),"Subtype"]=="TERT"], "SV", sep=";")
## We don't consider deletions in MYCN as drivers
matrix.early.late["MYCN",][matrix.early.late["MYCN",]==";DEL;Known_driver"] <- ""
## adjust heterogeneous cases
matrix.early.late["TERT","NBE57"] <- "SV"
matrix.early.late["MYCN","NBE57"] <- "AMP"
matrix.early.late["MYCN","NBE19"] <- "AMP"
matrix.early.late["TERT","NBE85"] <- "SV"

matrix.early.late["BRAF","NBE77"] <- "Subclonal;Known_driver"
matrix.early.late["MYCN",][matrix.early.late["MYCN",]==";DEL"] <- ""

## manually add as Known driver: ALK amplifications and F1174S, ATRX, MYCN, TERT mutations, except, deletions in CDKN2A, amplifications in CDK4 and MDM2
matrix.early.late["ALK",matrix.early.late["ALK",] %in% c(";AMP;SV", "AMP")] <- paste(matrix.early.late["ALK",matrix.early.late["ALK",] %in% c(";AMP;SV", "AMP", "Early or late clonal;")], "Known_driver", sep=";")
matrix.early.late["MYCN", matrix.early.late["MYCN",] %in% c(";AMP", "AMP")] <- paste(matrix.early.late["MYCN", matrix.early.late["MYCN",] %in% c(";AMP", "AMP")], "Known_driver", sep=";")
matrix.early.late["TERT", matrix.early.late["TERT",] %in% c("SV", ";SV")] <- paste(matrix.early.late["TERT", matrix.early.late["TERT",] %in% c(";SV", "SV")], "Known_driver", sep=";")
matrix.early.late["ATRX", matrix.early.late["ATRX",] %in% c(";DEL",  ";DUP")] <- paste(matrix.early.late["ATRX", matrix.early.late["ATRX",]  %in% c(";DEL" , ";DUP")], "Known_driver", sep=";")
matrix.early.late["CDKN2A", matrix.early.late["CDKN2A",] %in% c(";DEL", ";DEL;SV")] <- paste(matrix.early.late["CDKN2A", matrix.early.late["CDKN2A",] %in% c(";DEL", ";DEL;SV")], "Known_driver", sep=";")
matrix.early.late["CDK4", matrix.early.late["CDK4",] %in% c(";AMP", "AMP")] <- paste(matrix.early.late["CDK4", matrix.early.late["CDK4",] %in% c(";AMP", "AMP")], "Known_driver", sep=";")
matrix.early.late["MDM2", matrix.early.late["MDM2",] %in% c(";AMP", "AMP", ";AMP;SV")] <- paste(matrix.early.late["MDM2", matrix.early.late["MDM2",] %in% c(";AMP", "AMP", ";AMP;SV")], "Known_driver", sep=";")

## show whole chromosome info only if no partial gain
matrix.early.late["17", matrix.early.late["17q",]!=""] <- ""
matrix.early.late["1", matrix.early.late["1q",]!=""] <- ""
matrix.early.late["2", matrix.early.late["2p",]!=""] <- ""
matrix.early.late["7", matrix.early.late["7q",]!=""] <- ""

# ##############################################################################################################################################
# #### Count the number of hits/events
# 
# nr.hits <- apply(matrix.early.late, 2, function(x){
#   ## count chromosomal changes only as 2 events if there are at least 2 distinguishable events w.r.t. early late
#   tmp <- any(x[c("1", "1q", "7", "7q", "17", "17q", "2", "2p", "1p deletion", "11q deletion")] %in% c("early", "early tetraploid")) + 
#     any(x[c("1", "1q", "2", "2p", "7", "7q",  "17", "17q", "1p deletion", "11q deletion")] %in% c(" late", " late tetraploid"))
#   if(tmp==0){
#     tmp <- any(x[c("1", "1q", "7", "7q", "17", "17q", "2", "2p", "1p deletion", "11q deletion")] %in% c("early late", "early tetraploid late tetraploid"))
#   }
#   
#   tmp <- tmp + sum(sapply(x[setdiff(names(x), c("1", "1q", "7", "7q", "17", "17q", "2", "2p", "1p deletion", "11q deletion"))], function(y){
#     y <- strsplit(y, split=";")[[1]]
#     if(length(y)==0){
#       return(0)
#     }
#     y <- y[y!="Subclonal" & y!="Known_driver" & y!="NA"]
#     if(any(y!="")){
#       1
#     }else{
#       0
#     }
#   }))
#   tmp
# })
# 
# 
# 

