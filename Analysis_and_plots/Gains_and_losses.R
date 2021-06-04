##############################################################################################################################################
### Analyze gains and losses within the discovery cohort
##############################################################################################################################################

## homozygous-deletions due to artefacts need to be filtered
homdel.artifacts <- read.delim(paste0(meta.data, "artifact.homoDels.potentialArtifacts.txt"), sep="\t", header = T)

##############################################################################################################################################
## collect for each genomic position the number of tumors harboring a gain or loss at this position
gains <- GRanges()
losses <- gains

for(i in tumors.80x){
  ## ACEseq results files are named ...comb_pro...
  aceseq <- list.files(paste0(data.directory.80x,  i, "/", cnv.directory, "/"), pattern="comb_pro_extra")[1]
  if(is.na(aceseq)){next}
  
  purity. <- Extract.purity.ploidy.from.ACEseq(aceseq)$purity
  ploidy. <- Extract.purity.ploidy.from.ACEseq(aceseq)$ploidy
  
  ploidy. <- as.numeric(ploidy.)
  
  aceseq <-  read.delim(paste0(data.directory.80x,  i, "/", cnv.directory, "/", aceseq),  sep="\t", stringsAsFactors = F)
  
  aceseq$X.chromosome <- replace(aceseq$X.chromosome, aceseq$X.chromosome=="X", 23)
  aceseq$X.chromosome <- replace(aceseq$X.chromosome, aceseq$X.chromosome=="Y", 24)
  aceseq$X.chromosome <- as.numeric(aceseq$X.chromosome)
  aceseq <- aceseq[order(aceseq$X.chromosome),]
  
  
  ## filter artefacts
  for(j in 1:nrow(homdel.artifacts[1:12,])){
    tmp <- which(aceseq$X.chromosome==homdel.artifacts[j,1] & aceseq$start >= homdel.artifacts[j,2] & aceseq$end <= homdel.artifacts[j,3])
    if(length(tmp)>0){
      aceseq <- aceseq[-tmp,,drop=F]
    }
  }
  
  ## store the copy number info
  aceseq <-  GRanges(seqnames=aceseq[,1],
                     ranges = IRanges(aceseq[,2], aceseq[,3]),
                     CNV = aceseq[,13],
                     tumor = i,
                     gain.loss=aceseq[,ncol(aceseq)],
                     cov.ratio = aceseq[,"tcnMeanRaw"])
  
  
  gains <- c(gains, aceseq[(aceseq$gain.loss %in% c("gain", "LOHgain", "DUP", "DUP;LOH", "TCNneutral;DUP")| ((!is.na(aceseq$CNV) & as.numeric(aceseq$CNV) > (ploidy.)) & aceseq$cov.ratio>1.1)) |
                             (!is.na(aceseq$cov.ratio) && !is.na(aceseq$gain.loss)&& (aceseq$gain.loss == "sub" & aceseq$cov.ratio>1)),])
  
  losses <- c(losses,  aceseq[(aceseq$gain.loss %in% c("homozygousDel", "loss", "LOH", "cnLOH", "DEL", "HomoDel", "DEL;LOH") | (!is.na(aceseq$CNV) & as.numeric(aceseq$CNV) < (ploidy.) & aceseq$cov.ratio<.9)) |
                                (!is.na(aceseq$cov.ratio) && !is.na(aceseq$gain.loss)&& ( aceseq$gain.loss == "sub" & aceseq$cov.ratio<1)),])
}


gains$y<- 1
gains$Type <- "Gain"
losses$y <- -1
losses$Type <- "Loss"

data(ideoCyto, package = "biovizBase") ## data for ideogram
## use the lengths of the chromosomes from reference data
seqlengths(gains)[1:24] <- as.numeric(seqlengths(ideoCyto$hg19)[paste0("chr",c(1:22, "X", "Y"))])

setdiff(tumors.80x[telomere.classification.80x[tumors.80x]!="None"], unique(gains$tumor))
setdiff(names(telomere.classification.80x[telomere.classification.80x!="None"]), unique(losses$tumor))

## sum up the number of gains of individual tumors
gains <- coverage(gains)
start <- c()
end <- c()
## extract the positions (start end of each gain)
for(i in 1:length(gains)){
  position <- unlist(cumsum(unlist(runLength(gains[[i]]))))
  start <- c(start,1,position[-length(position)])
  end <- c(end,position)
}

## convert to GRanges object
gain.granges <- data.frame(Chr=rep(names(gains), sapply(gains, function(x) length(runLength(x)))), Strand = rep("*", length(unlist(runLength(gains)))),
                           Start = start, End = end,
                           Coverage = as.numeric(unlist(runValue(gains))))

gain.granges <- GRanges(seqnames = gain.granges$Chr, ranges = IRanges(start = gain.granges$Start, end = gain.granges$End), strand = 
                          gain.granges$Strand, Coverage = gain.granges$Coverage, Type = "gain")

gain.granges <- keepSeqlevels(gain.granges, as.character(1:24))
## use reference lengths of chromosomes
seqlengths(gain.granges) <- as.numeric(seqlengths(ideoCyto$hg19)[paste0("chr", c(1:22, "X", "Y"))])

## same for losses
losses <- coverage(losses)
start <- c()
end <- c()
for(i in 1:length(losses)){
  position <- unlist(cumsum(unlist(runLength(losses[[i]]))))
  if(length(position)==0){next}
  start <- c(start,1,position[-length(position)])
  end <- c(end,position)
}


loss.granges <- data.frame(Chr=rep(names(losses), sapply(losses, function(x) length(runLength(x)))), Strand = rep("*", length(unlist(runLength(losses)))),
                           Start = start, End = end,
                           Coverage = as.numeric(unlist(runValue(losses))))

loss.granges <- GRanges(seqnames = loss.granges$Chr, ranges = IRanges(start = loss.granges$Start, end = loss.granges$End), strand = 
                          loss.granges$Strand, Coverage = -loss.granges$Coverage, Type = "loss")

loss.granges <- keepSeqlevels(loss.granges, as.character(1:24)[as.character(1:24) %in% seqlevels(loss.granges)])

seqlengths(loss.granges)<- as.numeric(seqlengths(ideoCyto$hg19)[paste0("chr", c(1:22, "X", "Y"))[as.character(1:24) %in% seqlevels(loss.granges)]])


##############################################################################################################################################
## Quantify the number of chromosomes harboring a gain or loss, accounting for clonality state

### How many chromosomes do carry clonal or subclonal gains/losses?

chromosomes.with.clonal.gains <- c()
chromosomes.with.clonal.losses <- chromosomes.with.clonal.gains
chromosomes.with.subclonal.gains <- chromosomes.with.clonal.gains
chromosomes.with.subclonal.losses <- chromosomes.with.clonal.gains

for(i in tumors.80x){
  ## ACEseq results files are named ...comb_pro...
  aceseq <- list.files(paste0(data.directory.80x,  i, "/", cnv.directory, "/"), pattern="comb_pro_extra")[1]
  ploidy. <- strsplit(aceseq, split="extra")[[1]][2]
  ploidy. <- strsplit(ploidy., split="_")[[1]][1]
  ploidy. <- as.numeric(ploidy.)
  
  if(is.na(aceseq)){next}
  aceseq <-  read.delim(paste0(data.directory.80x,  i, "/", cnv.directory, "/", aceseq),  sep="\t", stringsAsFactors = F)
  
  aceseq$X.chromosome <- replace(aceseq$X.chromosome, aceseq$X.chromosome=="X", 23)
  aceseq$X.chromosome <- replace(aceseq$X.chromosome, aceseq$X.chromosome=="Y", 24)
  aceseq$X.chromosome <- as.numeric(aceseq$X.chromosome)
  aceseq <- aceseq[order(aceseq$X.chromosome),]
  
  ## filter artefacts
  for(j in 1:nrow(homdel.artifacts[1:12,])){
    tmp <- which(aceseq$X.chromosome==homdel.artifacts[j,1] & aceseq$start >= homdel.artifacts[j,2] & aceseq$end <= homdel.artifacts[j,3])
    if(length(tmp)>0){
      aceseq <- aceseq[-tmp,,drop=F]
    }
  }
  
  ## store the copy number info
  aceseq <-  GRanges(seqnames=aceseq[,1],
                     ranges = IRanges(aceseq[,2], aceseq[,3]),
                     CNV = aceseq[,"TCN"],
                     tumor = i,
                     gain.loss=aceseq[,ncol(aceseq)],
                     cov.ratio = aceseq[,"tcnMeanRaw"])
  
  ## only store gains and losses of a length of at least 10^6 bp
  aceseq <- aceseq[aceseq@ranges@width>=10^6,]
  aceseq <- aceseq[!is.na(aceseq$gain.loss),]
  aceseq.gains <- aceseq[!is.na(aceseq$gain.loss) & (aceseq$gain.loss %in% c("gain", "LOHgain", "DUP", "DUP;LOH", "TCNneutral;DUP") | (!is.na(as.numeric(aceseq$CNV)) & as.numeric(aceseq$CNV) > ploidy. & aceseq$cov.ratio>1.1)) & 
                           (!is.na(aceseq$CNV) & aceseq$CNV !="sub"),]
  chromosomes.with.clonal.gains <- c(chromosomes.with.clonal.gains, length(unique(seqnames(aceseq.gains))))
  names(chromosomes.with.clonal.gains)[length(chromosomes.with.clonal.gains)] <- i
  
  aceseq.subclonal.gains <- aceseq[(!is.na(aceseq$cov.ratio) & !is.na(aceseq$gain.loss)& 
                                      (!is.na(aceseq$CNV) & aceseq$CNV == "sub" & aceseq$cov.ratio>1.1)),]
  
  chromosomes.with.subclonal.gains <- c(chromosomes.with.subclonal.gains, length(unique(seqnames(aceseq.subclonal.gains))))
  
  aceseq.losses <- aceseq[(aceseq$gain.loss %in% c("homozygousDel", "loss", "LOH", "cnLOH", "DEL", "HomoDel", "DEL;LOH") | (aceseq$CNV!="sub" & as.numeric(aceseq$CNV) < ploidy. & aceseq$cov.ratio<0.9)) &
                            (!is.na(aceseq$CNV) & aceseq$CNV !="sub"),]
  
  chromosomes.with.clonal.losses <- c(chromosomes.with.clonal.losses,  length(unique(seqnames(aceseq.losses))))
  
  aceseq.subclonal.losses <- aceseq[(!is.na(aceseq$cov.ratio) & !is.na(aceseq$gain.loss)) &(!is.na(aceseq$CNV) &  
                                                                                              aceseq$CNV == "sub" & aceseq$cov.ratio<0.9),]
  
  chromosomes.with.subclonal.losses <- c(chromosomes.with.subclonal.losses,  length(unique(seqnames(aceseq.subclonal.losses))))
}

