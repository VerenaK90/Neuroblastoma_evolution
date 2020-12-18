##############################################################################################################################################
## Estimate the sampling bias of clonal mutations from tumor pairs

##############################################################################################################################################
## load libraries
library(bedr)

##############################################################################################################################################

## Store the clonal mutations and the false positive clonal mutations due to tissue sampling

## number of clonal mutations
clonal.mutations.all <- rep(0, length(primary.of.tumor.pairs))
names(clonal.mutations.all) <- primary.of.tumor.pairs
## number of false positive clonal mutations due to sampling
clonal.mutations.false.positives <- rep(0, length(primary.of.tumor.pairs))
names(clonal.mutations.false.positives) <- primary.of.tumor.pairs

## indicator vectors for the copy number and the B allele count

copy.number.indicator <- 1:6

## iterate through the tumors, then iterate through the different copy number states

## outmost loop: tumors
for(i in 1:length(primary.of.tumor.pairs)){
  

  prim <- primary.of.tumor.pairs[i]
  relapse <- relapse.of.tumor.pairs[i]
  
  ## ploidy/purity:
  aceseq.prim <- list.files(paste0(data.directory, "/", prim, "/ACEseq/"), pattern="comb_pro_extra")[1]
  if(is.na(aceseq.prim)){next}
  ploidy.prim <- sapply(aceseq.prim, function(x){strsplit(x, split="extra")[[1]][2]})
  purity.prim <- sapply(ploidy.prim, function(x){strsplit(x, split="_")[[1]][2]})
  purity.prim <- as.numeric(sapply(purity.prim, function(x){strsplit(x, split=".txt")[[1]][1]}))
  ploidy.prim <- as.numeric(sapply(ploidy.prim, function(x){strsplit(x, split="_")[[1]][1]}))
  
  aceseq.rec <- list.files(paste0(data.directory, "/", relapse, "/ACEseq/"), pattern="comb_pro_extra")[1]
  if(is.na(aceseq.rec)){next}
  ploidy.rec <- sapply(aceseq.rec, function(x){strsplit(x, split="extra")[[1]][2]})
  purity.rec <- sapply(ploidy.rec, function(x){strsplit(x, split="_")[[1]][2]})
  purity.rec <- as.numeric(sapply(purity.rec, function(x){strsplit(x, split=".txt")[[1]][1]}))
  ploidy.rec <- as.numeric(sapply(ploidy.rec, function(x){strsplit(x, split="_")[[1]][1]}))
  
  ## read in the mutation file
  files.prim <- list.files(paste0(data.directory, "/", prim, "/SNVs/"), pattern="somatic_snvs_conf_8_to_10")[1]
  
  mutations.prim <- read.vcf(paste0(data.directory, "/", prim, "/SNVs/", files.prim))
  readcounts.prim <- t(sapply(mutations.prim$vcf$INFO, function(x){
    x <- strsplit(x, split=";")[[1]][6]
    x <- strsplit(x, split="=")[[1]][2]
    x <- as.numeric(strsplit(x, split=",")[[1]])
    x <- unname(x)
    c(sum(x[c(1,2)]), sum(x[c(3,4)]))
  }))
  rownames(readcounts.prim) <- paste(mutations.prim$vcf$GENE, mutations.prim$vcf$POS, sep=".")
  
  ### the same for the relapse sample
  files.rec <- list.files(paste0(data.directory, "/", relapse, "/SNVs/"), pattern="somatic_snvs_conf_8_to_10")[1]
  mutations.rec <- read.vcf(paste0(data.directory, "/", relapse, "/SNVs/", files.prim))
  readcounts.rec <- t(sapply(mutations.rec$vcf$INFO, function(x){
    x <- strsplit(x, split=";")[[1]][6]
    x <- strsplit(x, split="=")[[1]][2]
    x <- as.numeric(strsplit(x, split=",")[[1]])
    x <- unname(x)
    c(sum(x[c(1,2)]), sum(x[c(3,4)]))
  }))
  rownames(readcounts.rec) <- paste(mutations.rec$vcf$GENE, mutations.rec$vcf$POS, sep=".")
  
  ## extract only the sites of the primary tumor, since we simply want to know how many of them are lost in the relapse
  
  readcounts.rec. <-  matrix(0, nrow=nrow(readcounts.prim), ncol=2,
                             dimnames = list(rownames(readcounts.prim)))
  
  readcounts.rec.[intersect(rownames(readcounts.rec), rownames(readcounts.prim)),] <- readcounts.rec[intersect(rownames(readcounts.rec), rownames(readcounts.prim)),]

  readcounts.rec <- readcounts.rec.
  rm(readcounts.rec.)
  
  
  copy.number.info.prim <- read.delim(file=paste0(data.directory, "/", prim, "/ACEseq/", aceseq.prim), sep="\t", stringsAsFactors = F)
  ## obtain the coverage ratios for the mutations of interest. These are the mutations present in the primary tumor
  coverage.ratios.prim <- apply(mutations.prim$vcf, 1, function(x){
    x <- unlist(x)
    tmp <- copy.number.info.prim[copy.number.info.prim$X.chromosome==x[1],]
    tmp <- tmp[as.numeric(tmp$start) <= as.numeric(x[2]) & as.numeric(tmp$end) >= as.numeric(x[2]),]
    if(nrow(tmp)>0){
      unlist(tmp[28])
    }else{
      NA
    }
  })
  ## and the B-allele frequencies
  bafs.prim <- apply(mutations.prim$vcf, 1, function(x){
    x <- unlist(x)
    tmp <- copy.number.info.prim[copy.number.info.prim$X.chromosome==x[1],]
    tmp <- tmp[as.numeric(tmp$start) <= as.numeric(x[2]) & as.numeric(tmp$end) >= as.numeric(x[2]),]
    if(nrow(tmp)>0){
      unlist(tmp[30])
    }else{
      NA
    }
  })
  ## and the genotype
  genotype.prim <- apply(mutations.prim$vcf, 1, function(x){
    x <- unlist(x)
    tmp <- copy.number.info.prim[copy.number.info.prim$X.chromosome==x[1],]
    tmp <- tmp[as.numeric(tmp$start) <= as.numeric(x[2]) & as.numeric(tmp$end) >= as.numeric(x[2]),]
    if(nrow(tmp)>0){
      unlist(tmp[36])
    }else{
      NA
    }
  })
  ## and the suggested copy number
  tcn.prim <- apply(mutations.prim$vcf, 1, function(x){
    x <- unlist(x)
    tmp <- copy.number.info.prim[copy.number.info.prim$X.chromosome==x[1],]
    tmp <- tmp[as.numeric(tmp$start) <= as.numeric(x[2]) & as.numeric(tmp$end) >= as.numeric(x[2]),]
    if(nrow(tmp)>0){
      as.numeric(unlist(tmp[37]))
    }else{
      NA
    }
  })
  
  depth.prim <- median(rowSums(readcounts.prim), na.rm=T)
  
   depth.rec <- median(rowSums(readcounts.rec), na.rm=T)
  

    #######################################################################
    ## second level of iteration: go through all copy number states. 
    
    for(k in unique(copy.number.indicator)){
      
        
        ## determine the positions of the clonal peaks
        haploid.prob.clonal.prim <- purity.prim/(purity.prim*k + (1-purity.prim)*2)
        
        haploid.prob.clonal.rec <- purity.rec/(purity.rec*k + (1-purity.rec)*2)
        
       
        ## tolerance +/- 0.1 from expected coverage ratio at triploid fraction
        expected.coverage.ratio.prim <- (k*purity.prim + (1-purity.prim)*2)/(ploidy.prim*purity.prim+(1-purity.prim)*2)
        expected.coverage.ratio.rec <- (k*purity.rec + (1-purity.rec)*2)/(ploidy.rec*purity.rec+(1-purity.rec)*2)
        
        ## how many mutations lie on the respective allele?
        ## 1st: choose all mutations that lie on the specific aneuploidy loci in the primary tumor
        mutations.at.copy.number.prim <-  readcounts.prim[((coverage.ratios.prim>(k/expected.coverage.ratio.prim-0.1) & 
                                                              coverage.ratios.prim < (expected.coverage.ratio.prim+0.1) & 
                                                              !is.na(coverage.ratios.prim)) |
                                                   tcn.prim==k & !is.na(tcn.prim)),,drop=F]
        
        
        mutations.at.copy.number.rec <-  readcounts.rec[((coverage.ratios.prim>(expected.coverage.ratio.prim-0.1) & 
                                                            coverage.ratios.prim < (expected.coverage.ratio.prim+0.1) & 
                                                              !is.na(coverage.ratios.prim)) |
                                                             tcn.prim==k & !is.na(tcn.prim))  ,,drop=F]
        

        
        if(nrow(mutations.at.copy.number.prim)==0){next}
        
        ## note: deal with haploid and heterozygous diploid cases separately
        
        ## 2nd we now count all mutations that lie on the right hand side of the haploid clonal peak 

        mutations.at.copy.number.prim <- mutations.at.copy.number.prim[mutations.at.copy.number.prim[,2]/rowSums(mutations.at.copy.number.prim)>=haploid.prob.clonal.prim,,drop=F]
        mutations.at.copy.number.prim <- mutations.at.copy.number.prim[!is.na(rowSums(mutations.at.copy.number.prim)),,drop=F]
        
        ## just check whether these mutations disappeared as a crude measure of false positives
        mutations.at.copy.number.rec <- mutations.at.copy.number.rec[!is.na(rowSums(mutations.at.copy.number.rec)),,drop=F]
        
        mutations.at.copy.number.rec <- mutations.at.copy.number.rec[intersect(rownames(mutations.at.copy.number.prim),
                                                                               rownames(mutations.at.copy.number.rec)),,drop=F]
        
        mutations.at.copy.number.prim <- mutations.at.copy.number.prim[rowSums(mutations.at.copy.number.prim)>0,,drop=F]
        mutations.at.copy.number.rec <- mutations.at.copy.number.rec[rowSums(mutations.at.copy.number.rec)>0,,drop=F]
        
        if(nrow(mutations.at.copy.number.prim)==0){next}
        
        clonal.mutations.all[i] <- clonal.mutations.all[i] + nrow(mutations.at.copy.number.prim)
        clonal.mutations.false.positives[i] <- clonal.mutations.false.positives[i] + nrow(mutations.at.copy.number.prim) - nrow(mutations.at.copy.number.rec)
   
    }
  
}
