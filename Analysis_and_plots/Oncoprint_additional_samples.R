
##############################################################################################################################################

## in addition, add chromosomal gains at early and late stages
tumors.with.17q <- rownames(sample.information.30x)[sample.information.30x$`17q.gain`=="T"]
tumors.with.whole.17.gain <- rownames(sample.information.30x)[sample.information.30x$`17.gain`=="T" & sample.information.30x$`17q.gain`=="F"]
tumors.with.7q <- rownames(sample.information.30x)[sample.information.30x$`7q.gain`=="T"]
tumors.with.whole.7.gain <- rownames(sample.information.30x)[sample.information.30x$`7.gain`=="T" & sample.information.30x$`7q.gain`=="F"]
tumors.with.1p.loss <- rownames(sample.information.30x)[sample.information.30x$`1p.deletion`=="T"]
tumors.with.1q.gain <- rownames(sample.information.30x)[sample.information.30x$`1q.gain`=="T"]
tumors.with.whole.1.gain <- rownames(sample.information.30x)[sample.information.30x$`1.gain`=="T" & sample.information.30x$`1q.gain`=="F"]
tumors.with.2p <- rownames(sample.information.30x)[sample.information.30x$`2p.gain`=="T"]
tumors.with.whole.2.gain <- rownames(sample.information.30x)[sample.information.30x$`2.gain`=="T" & sample.information.30x$`2p.gain`=="F"]
tumors.with.11q.loss <- rownames(sample.information.30x)[sample.information.30x$`11q.deletion`=="T"]

not.dateable <- data.frame(Gain = "17q",
                           Sample = rownames(sample.information.30x)[sample.information.30x$`17q.gain`=="T" & sample.information.30x$`17q.gain.<.5`=="F"])
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "17",
                                 Sample = rownames(sample.information.30x)[sample.information.30x$`17.gain`=="T" & sample.information.30x$`17.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "7q",
                                 Sample = rownames(sample.information.30x)[sample.information.30x$`7q.gain`=="T" & sample.information.30x$`7q.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "7",
                                 Sample = rownames(sample.information.30x)[sample.information.30x$`7.gain`=="T" & sample.information.30x$`7.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "1q",
                                 Sample = rownames(sample.information.30x)[sample.information.30x$`1q.gain`=="T" & sample.information.30x$`1q.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "1",
                                 Sample = rownames(sample.information.30x)[sample.information.30x$`1.gain`=="T" & sample.information.30x$`1.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "2p",
                                 Sample = rownames(sample.information.30x)[sample.information.30x$`2p.gain`=="T" & sample.information.30x$`2p.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "2",
                                 Sample = rownames(sample.information.30x)[sample.information.30x$`2.gain`=="T" & sample.information.30x$`2.gain.<.5`=="F"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "11q",
                                 Sample = rownames(sample.information.30x)[sample.information.30x$`11q.deletion`=="T"]))
not.dateable <- rbind(not.dateable,
                      data.frame(Gain = "1p",
                                 Sample = rownames(sample.information.30x)[sample.information.30x$`1p.deletion`=="T"]))



## to this end, read in the timing information
load(paste0(rdata.directory, "MRCA_timing.RData"))
load(paste0(rdata.directory,"Clonal_mutations_different_ploidies.RData"))

subclonal.gains <- data.frame(CHROM=c("1q", "17q", "17", "17q", "7", "1q"), SAMPLE=c("B087wgs_23206", "B087koeln_23229", "B087wgs_24719", "K09R-SPVTXJ",
                                                                         "B087wgs_25262", "K09R-KJ62ZX"))


## read in annotation
gains.losses.eca <- data.frame(CHROM=c(), SAMPLE=c())
gains.losses.eca <- rbind(gains.losses.eca, 
                          data.frame(CHROM=rep("1", sum(sample.information.30x$Timing_whole_1 %in% c("ECA", "ECA/MRCA"))), 
                                     SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_whole_1 %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, 
                          data.frame(CHROM="1p", SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_1p %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="1q", 
                                                       SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_1q %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="2", 
                                                       SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_whole_2 %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="2p", 
                                                       SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_2p %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="7", 
                                                       SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_whole_7 %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="7q", 
                                                       SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_7q %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="11q", 
                                                       SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_11q %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="17", 
                                                       SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_whole_17 %in% c("ECA", "ECA/MRCA"),])))
gains.losses.eca <- rbind(gains.losses.eca, data.frame(CHROM="17q", 
                                                       SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_17q %in% c("ECA", "ECA/MRCA"),])))

gains.losses.early.tetraploid <- data.frame(CHROM=c(), SAMPLE=c())
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, 
                                       data.frame(CHROM=rep(1, sum(sample.information.30x$Timing_whole_1 %in% c("early tetraploid"))), SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_whole_1 %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, 
                                       data.frame(CHROM=rep("1p", sum(sample.information.30x$Timing_1p %in% c("early tetraploid"))), SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_1p %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("1q", sum(sample.information.30x$Timing_1q %in% c("early tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_1q %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("2", sum(sample.information.30x$Timing_whole_2 %in% c("early tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_whole_2 %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("2p", sum(sample.information.30x$Timing_2p %in% c("early tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_2p %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("7", sum(sample.information.30x$Timing_whole_7 %in% c("early tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_whole_7 %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("7q", sum(sample.information.30x$Timing_7q %in% c("early tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_7q %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, 
                                       data.frame(CHROM=rep("11q", sum(sample.information.30x$Timing_11q %in% c("early tetraploid"))), SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_11q %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("17", sum(sample.information.30x$Timing_whole_17 %in% c("early tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_whole_17 %in% c("early tetraploid"),])))
gains.losses.early.tetraploid <- rbind(gains.losses.early.tetraploid, data.frame(CHROM=rep("17q", sum(sample.information.30x$Timing_17q %in% c("early tetraploid"))), 
                                                                                 SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_17q %in% c("early tetraploid"),])))


gains.losses.mrca <- data.frame(CHROM=1, SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_whole_1 %in% c("MRCA", "ECA/MRCA"),]))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="1p", 
                                                         SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_1p %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="1q", 
                                                         SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_1q %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="2", 
                                                         SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_whole_2 %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="2p", 
                                                         SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_2p %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="7", 
                                                         SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_whole_7 %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="7q", 
                                                         SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_7q %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="11q", 
                                                         SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_11q %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="17", 
                                                         SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_whole_17 %in% c("MRCA", "ECA/MRCA"),])))
gains.losses.mrca <- rbind(gains.losses.mrca, data.frame(CHROM="17q", 
                                                         SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_17q %in% c("MRCA", "ECA/MRCA"),])))


gains.losses.late.tetraploid <- data.frame(CHROM=c(), SAMPLE=c())
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, 
                                      data.frame(CHROM=rep(1, sum(sample.information.30x$Timing_whole_1 %in% c("late tetraploid"))), SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_whole_1 %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("1q", sum(sample.information.30x$Timing_1q %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_1q %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("1p", sum(sample.information.30x$Timing_1p %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_1p %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("2", sum(sample.information.30x$Timing_whole_2 %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_whole_2 %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("2p", sum(sample.information.30x$Timing_2p %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_2p %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("7", sum(sample.information.30x$Timing_whole_7 %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_whole_7 %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("7q", sum(sample.information.30x$Timing_7q %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_11q %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("11q", sum(sample.information.30x$Timing_11q %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_7q %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("17", sum(sample.information.30x$Timing_whole_17 %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_whole_17 %in% c("late tetraploid"),])))
gains.losses.late.tetraploid <- rbind(gains.losses.late.tetraploid, data.frame(CHROM=rep("17q", sum(sample.information.30x$Timing_17q %in% c("late tetraploid"))), 
                                                                               SAMPLE=rownames(sample.information.30x[sample.information.30x$Timing_17q %in% c("late tetraploid"),])))


##############################################################################################################################################
#### design a mutation matrix, where the early/late information is included
matrix.early.late.30x <- matrix("", nrow=10,
                            ncol=length(tumors.30x), dimnames=list(c("7", "2", "1", "17", "7q", "2p", "1q", "17q", "11q", "1p"),
                                                               tumors.30x))


for(i in rownames(matrix.early.late.30x)){
  matrix.early.late.30x[i,as.character(gains.losses.eca[gains.losses.eca$CHROM==i,]$SAMPLE)] <- "early"
  matrix.early.late.30x[i,as.character(gains.losses.mrca[gains.losses.mrca$CHROM==i,]$SAMPLE)] <- paste(matrix.early.late.30x[i,as.character(gains.losses.mrca[gains.losses.mrca$CHROM==i,]$SAMPLE)], "late", sep=" ")
  matrix.early.late.30x[i,as.character(subclonal.gains[subclonal.gains$CHROM==i,]$SAMPLE)] <- "Subclonal"
  matrix.early.late.30x[i,as.character(gains.losses.early.tetraploid[gains.losses.early.tetraploid$CHROM==i,]$SAMPLE)] <- "early tetraploid"
  matrix.early.late.30x[i,as.character(gains.losses.late.tetraploid[gains.losses.late.tetraploid$CHROM==i,]$SAMPLE)] <- paste(matrix.early.late.30x[i,as.character(gains.losses.late.tetraploid[gains.losses.late.tetraploid$CHROM==i,]$SAMPLE)], "late tetraploid", sep=" ")
    
}


matrix.early.late.30x[matrix.early.late.30x==0] <- ""

## show whole chromosome info only if no partial gain
matrix.early.late.30x["17", matrix.early.late.30x["17q",]!=""] <- ""
matrix.early.late.30x["1", matrix.early.late.30x["1q",]!=""] <- ""
matrix.early.late.30x["2", matrix.early.late.30x["2p",]!=""] <- ""
matrix.early.late.30x["7", matrix.early.late.30x["7q",]!=""] <- ""
