##############################################################################################################################################
## For 2 tumors, we have sequencing data from the primary and the relapse. We use them to estimate the sampling bias in the call of clonal mutations.

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
  aceseq.prim <- list.files(paste0(data.directory.discovery, "/", prim, "/", cnv.directory, "/"), pattern="comb_pro_extra")[1]
  if(is.na(aceseq.prim)){next}
  ploidy.prim <- Extract.purity.ploidy.from.ACEseq(aceseq.prim)$ploidy
  purity.prim <- Extract.purity.ploidy.from.ACEseq(aceseq.prim)$purity

  aceseq.rec <- list.files(paste0(data.directory.discovery, "/", relapse, "/", cnv.directory, "/"), pattern="comb_pro_extra")[1]
  if(is.na(aceseq.rec)){next}
  ploidy.rec <- Extract.purity.ploidy.from.ACEseq(aceseq.rec)$ploidy
  purity.rec <- Extract.purity.ploidy.from.ACEseq(aceseq.rec)$purity

  
  ## read in the mutation file
  files.prim <- list.files(paste0(data.directory.discovery, "/", prim, "/", snv.directory, "/"), pattern="somatic_snvs_conf_8_to_10")[1]
  
  mutations.prim <- read.vcf(paste0(data.directory.discovery, "/", prim, "/", files.prim))
  readcounts.prim <- Extract.info.from.vcf(mutations.prim, info="readcounts", type="snvs", mutationcaller = "DKFZ")

  ### the same for the relapse sample
  files.rec <- list.files(paste0(data.directory.discovery, "/", relapse,  "/", snv.directory, "/"), pattern="somatic_snvs_conf_8_to_10")[1]
  mutations.rec <- read.vcf(paste0(data.directory.discovery, "/", relapse,  "/", snv.directory, "/", files.prim))
  readcounts.rec <-  Extract.info.from.vcf(mutations.rec, info="readcounts", type="snvs", mutationcaller = "DKFZ")
  
  ## extract only the sites of the primary tumor, since we simply want to know how many of them are lost in the relapse
  
  readcounts.rec. <-  matrix(0, nrow=nrow(readcounts.prim), ncol=2,
                             dimnames = list(rownames(readcounts.prim)))
  
  readcounts.rec.[intersect(rownames(readcounts.rec), rownames(readcounts.prim)),] <- readcounts.rec[intersect(rownames(readcounts.rec), rownames(readcounts.prim)),]
  
  readcounts.rec <- readcounts.rec.
  rm(readcounts.rec.)
  
  
  copy.number.info.prim <- read.delim(file=paste0(data.directory.discovery, "/", prim, "/", cnv.directory, "/", aceseq.prim), sep="\t", stringsAsFactors = F)
  ## obtain the coverage ratios/BAF/genotype/TCN for the mutations of interest. These are the mutations present in the primary tumor
  cnv.info.per.mutation.prim <- Extract.copy.number.info.per.SSNV(mutations.prim, copy.number.info.prim)
  ## obtain the coverage ratios at mutated loci
  coverage.ratios.prim <- cnv.info.per.mutation.prim$coverage.ratio
  bafs.prim <- cnv.info.per.mutation.prim$baf
  genotype.prim <- cnv.info.per.mutation.prim$genotype
  tcn.prim <- cnv.info.per.mutation.prim$tcn
  
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

## if re-running make sure that manually estimated purity values are okay
save(clonal.mutations.all, clonal.mutations.false.positives,
     file=paste0(rdata.directory, "Assess_false_positives_due_to_sampling.RData"))