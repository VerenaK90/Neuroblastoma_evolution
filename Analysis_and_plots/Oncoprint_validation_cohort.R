##############################################################################################################################################
## read in the driver mutation information
source("Settings.R")

##############################################################################################################################################

## in addition, add chromosomal gains at early and late stages
tumors.with.17q <- rownames(sample.information.validation)[sample.information.validation$`17q.gain`=="T"]
tumors.with.whole.17.gain <- rownames(sample.information.validation)[sample.information.validation$`17.gain`=="T" & sample.information.validation$`17q.gain`=="F"]
tumors.with.7q <- rownames(sample.information.validation)[sample.information.validation$`7q.gain`=="T"]
tumors.with.whole.7.gain <- rownames(sample.information.validation)[sample.information.validation$`7.gain`=="T" & sample.information.validation$`7q.gain`=="F"]
tumors.with.1p.loss <- rownames(sample.information.validation)[sample.information.validation$`1p.deletion`=="T"]
tumors.with.1q.gain <- rownames(sample.information.validation)[sample.information.validation$`1q.gain`=="T"]
tumors.with.whole.1.gain <- rownames(sample.information.validation)[sample.information.validation$`1.gain`=="T" & sample.information.validation$`1q.gain`=="F"]
tumors.with.2p <- rownames(sample.information.validation)[sample.information.validation$`2p.gain`=="T"]
tumors.with.whole.2.gain <- rownames(sample.information.validation)[sample.information.validation$`2.gain`=="T" & sample.information.validation$`2p.gain`=="F"]
tumors.with.11q.loss <- rownames(sample.information.validation)[sample.information.validation$`11q.deletion`=="T"]

not.dateable <- data.frame(Gain = "17q",
                           Sample = rownames(sample.information.validation)[sample.information.validation$`17q.gain`=="T" & sample.information.validation$`17q.gain.<.5`=="F"])
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "17",
                                 Sample = rownames(sample.information.validation)[sample.information.validation$`17.gain`=="T" & sample.information.validation$`17.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "7q",
                                 Sample = rownames(sample.information.validation)[sample.information.validation$`7q.gain`=="T" & sample.information.validation$`7q.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "7",
                                 Sample = rownames(sample.information.validation)[sample.information.validation$`7.gain`=="T" & sample.information.validation$`7.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "1q",
                                 Sample = rownames(sample.information.validation)[sample.information.validation$`1q.gain`=="T" & sample.information.validation$`1q.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "1",
                                 Sample = rownames(sample.information.validation)[sample.information.validation$`1.gain`=="T" & sample.information.validation$`1.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "2p",
                                 Sample = rownames(sample.information.validation)[sample.information.validation$`2p.gain`=="T" & sample.information.validation$`2p.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "2",
                                 Sample = rownames(sample.information.validation)[sample.information.validation$`2.gain`=="T" & sample.information.validation$`2.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "11q",
                                 Sample = rownames(sample.information.validation)[sample.information.validation$`11q.deletion`=="T"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "1p",
                                 Sample = rownames(sample.information.validation)[sample.information.validation$`1p.deletion`=="T"]))


## manual annotation of subclonal gains

subclonal.gains <- data.frame(CHROM=c("1q", "17q", "17", "17q", "7", "1q"), SAMPLE=c("NBE114", "NBE115", "NBE124", "NBE127",
                                                                         "NBE129", "NBE135"))


## read in annotation
gains.losses.eca <- data.frame(CHROM=c(), SAMPLE=c())
gains.losses.eca <- rbind(gains.losses.eca, 
                          data.frame(CHROM=rep("1", sum(sample.information.validation$Timing_whole_1 %in% c("ECA", "ECA/MRCA"))), 
                                     SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_whole_1 %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, 
                          data.frame(CHROM="1p", SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_1p %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="1q", 
                                                       SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_1q %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="2", 
                                                       SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_whole_2 %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="2p", 
                                                       SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_2p %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="7", 
                                                       SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_whole_7 %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="7q", 
                                                       SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_7q %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="11q", 
                                                       SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_11q %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="17", 
                                                       SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_whole_17 %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="17q", 
                                                       SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_17q %in% c("ECA", "ECA/MRCA"),])))

gains.losses.early.tetraploid <- data.frame(CHROM=c(), SAMPLE=c())
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, 
                                       data.frame(CHROM=rep(1, sum(sample.information.validation$Timing_whole_1 %in% c("early tetraploid"))), SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_whole_1 %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, 
                                       data.frame(CHROM=rep("1p", sum(sample.information.validation$Timing_1p %in% c("early tetraploid"))), SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_1p %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("1q", sum(sample.information.validation$Timing_1q %in% c("early tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_1q %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("2", sum(sample.information.validation$Timing_whole_2 %in% c("early tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_whole_2 %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("2p", sum(sample.information.validation$Timing_2p %in% c("early tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_2p %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("7", sum(sample.information.validation$Timing_whole_7 %in% c("early tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_whole_7 %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("7q", sum(sample.information.validation$Timing_7q %in% c("early tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_7q %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, 
                                       data.frame(CHROM=rep("11q", sum(sample.information.validation$Timing_11q %in% c("early tetraploid"))), SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_11q %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("17", sum(sample.information.validation$Timing_whole_17 %in% c("early tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_whole_17 %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("17q", sum(sample.information.validation$Timing_17q %in% c("early tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_17q %in% c("early tetraploid"),])))


gains.losses.mrca <- data.frame(CHROM=1, SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_whole_1 %in% c("MRCA", "ECA/MRCA"),]))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="1p", 
                                                         SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_1p %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="1q", 
                                                         SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_1q %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="2", 
                                                         SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_whole_2 %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="2p", 
                                                         SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_2p %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="7", 
                                                         SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_whole_7 %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="7q", 
                                                         SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_7q %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="11q", 
                                                         SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_11q %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="17", 
                                                         SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_whole_17 %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="17q", 
                                                         SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_17q %in% c("MRCA", "ECA/MRCA"),])))


gains.losses.late.tetraploid <- data.frame(CHROM=c(), SAMPLE=c())
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, 
                                      data.frame(CHROM=rep(1, sum(sample.information.validation$Timing_whole_1 %in% c("late tetraploid"))), SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_whole_1 %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("1q", sum(sample.information.validation$Timing_1q %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_1q %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("1p", sum(sample.information.validation$Timing_1p %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_1p %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("2", sum(sample.information.validation$Timing_whole_2 %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_whole_2 %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("2p", sum(sample.information.validation$Timing_2p %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_2p %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("7", sum(sample.information.validation$Timing_whole_7 %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_whole_7 %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("7q", sum(sample.information.validation$Timing_7q %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_11q %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("11q", sum(sample.information.validation$Timing_11q %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_7q %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("17", sum(sample.information.validation$Timing_whole_17 %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_whole_17 %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("17q", sum(sample.information.validation$Timing_17q %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.validation[sample.information.validation$Timing_17q %in% c("late tetraploid"),])))


##############################################################################################################################################
#### design a mutation matrix, where the early/late information is included
matrix.early.late.validation <- matrix("", nrow=10,
                            ncol=length(tumors.validation), dimnames=list(c("7", "2", "1", "17", "7q", "2p", "1q", "17q", "11q", "1p"),
                                                               tumors.validation))


for(i in rownames(matrix.early.late.validation)){
  matrix.early.late.validation[i,as.character(gains.losses.eca[gains.losses.eca$CHROM==i,]$SAMPLE)] <- "early"
  matrix.early.late.validation[i,as.character(gains.losses.mrca[gains.losses.mrca$CHROM==i,]$SAMPLE)] <- paste(matrix.early.late.validation[i,as.character(gains.losses.mrca[gains.losses.mrca$CHROM==i,]$SAMPLE)], "late", sep=" ")
  matrix.early.late.validation[i,as.character(subclonal.gains[subclonal.gains$CHROM==i,]$SAMPLE)] <- "Subclonal"
  matrix.early.late.validation[i,as.character(gains.losses.early.tetraploid[gains.losses.early.tetraploid$CHROM==i,]$SAMPLE)] <- "early tetraploid"
  matrix.early.late.validation[i,as.character(gains.losses.late.tetraploid[gains.losses.late.tetraploid$CHROM==i,]$SAMPLE)] <- paste(matrix.early.late.validation[i,as.character(gains.losses.late.tetraploid[gains.losses.late.tetraploid$CHROM==i,]$SAMPLE)], "late tetraploid", sep=" ")
    
}


matrix.early.late.validation[matrix.early.late.validation==0] <- ""

## show whole chromosome info only if no partial gain
matrix.early.late.validation["17", matrix.early.late.validation["17q",]!=""] <- ""
matrix.early.late.validation["1", matrix.early.late.validation["1q",]!=""] <- ""
matrix.early.late.validation["2", matrix.early.late.validation["2p",]!=""] <- ""
matrix.early.late.validation["7", matrix.early.late.validation["7q",]!=""] <- ""
