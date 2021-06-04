##############################################################################################################################################
## Quantify the densitiy of amplified and non-amplified clonal mutations
##############################################################################################################################################
## specify sample pairs and source sampling bias file

##  primary tumor of tumor pairs
primary.of.tumor.pairs <- c("NBE11", "NBE51")
## and relapse tumor of tumor pairs
relapse.of.tumor.pairs <- c("NBE66", "NBE78")

source(paste0(custom.script.directory, "Estimate_sampling_bias_from_tumor_pairs.R"))
load(paste0(rdata.directory, "Assess_false_positives_due_to_sampling.RData"))
##############################################################################################################################################
## For each of the tumors, we want to know the number of clonal mutations that lie on one or multiple copies of a gained chromosome. 
## Clonal mutations on several copies were likely acquired prior to tumor initiation, given that the gain is clonal. 
## Since the number is only informative if corrected for the respective genome fraction, we need to report this as well.

## Store ACEseq estimates of purity and ploidy for each tumor
purity <- c()
ploidy <- purity

## Store the clonal mutations present on any number of alleles and store this information for each chromosome separately in a list.
## The roman numbers correspond to the number of alleles.

clonal.mutation.matrix. <- matrix(0, nrow=1 + 2 + 3 + 4 , ncol=length(c(tumors.80x, tumors.30x)),
                                 dimnames=list(c("monosomic", "disomic_I", "disomic_II", "trisomic_I", "trisomic_II", "trisomic_III",
                                                 "tetrasomic_I", "tetrasomic_II", "tetrasomic_III", "tetrasomic_IV"),
                                               c(tumors.80x, tumors.30x)))
## one matrix per autosome
clonal.mutation.matrix <- list()

for(i in 1:22){
  clonal.mutation.matrix[[i]] <- clonal.mutation.matrix.
  names(clonal.mutation.matrix)[i] <- paste0("chr",i)
}
rm(clonal.mutation.matrix.)


## analogously, we need a genome length table

segment.length.matrix <- clonal.mutation.matrix


## Iterate through the tumors, iterate through the autosomes and iterate through the different B allele frequencies

## First loop: tumors
for(i in c(tumors.80x, tumors.30x)){
  
  print(i)
  
  if(i %in% tumors.80x){
    data.directory <- data.directory.80x
  }else{
    data.directory <- data.directory.30x
  }
  
  ## Find ACEseq file:
  aceseq <- list.files(paste0(data.directory, "/", i, "/", cnv.directory, "/"), pattern="comb_pro_extra")[1]
  if(is.na(aceseq)){next}
  aceseq <- paste0(data.directory, "/", i, "/", cnv.directory, "/", aceseq)
  
  ## Find mutation file
  mutations <- list.files(paste0(data.directory, "/", i, "/", snv.directory, "/"), pattern="somatic_snvs_conf_8_to_10")[1]
  mutations <- paste0(data.directory, "/", i, "/", snv.directory, "/", mutations)
  
  ## Correct purity estimate for tumor NBE40 manually
  if(i == "NBE40"){
    purity.refit <- T
  }else{
    purity.refit <- F
  }
  clonal.mutations <- count.clonal.mutations(aceseq, mutations, purity.refit=purity.refit, chromosomes = c(1:22))
  
  purity[i] <- clonal.mutations$purity
  ploidy[i] <- clonal.mutations$ploidy
  
  for(j in 1:22){
    clonal.mutation.matrix[[j]][,i] <- clonal.mutations$clonal.mutation.matrix[,j]
    segment.length.matrix[[j]][,i] <- clonal.mutations$segment.length.matrix[,j]
  }
  
  copy.number.indicator <- clonal.mutations$copy.number.indicator
  B.allele.indicator <- clonal.mutations$B.allele.indicator

}

save(ploidy, purity, copy.number.indicator, B.allele.indicator, segment.length.matrix, clonal.mutation.matrix, 
     file=paste0(rdata.directory, "Clonal_mutations_different_ploidies.RData"))

load(paste0(rdata.directory, "Clonal_mutations_different_ploidies.RData"))


##############################################################################################################################################
## Compute the mutational density at each genomic fragment 

## The monosomic peak needs to be multiplied with two, since we only estimated the RHS. 

## We model mutation accumulation as a Poisson process. The number of mutations acquired on a piece of DNA depends on the per-base mutation rate (at this time) and the length of the piece
## using a chisquare-approximation, we obtain confidence intervals for the mutation time. I store these in a list

# and store the mutation time of the 1st event

mutation.time.most.likely <- list()
mutation.time.lower <- list()
mutation.time.upper <- list()

for(i in c(tumors.80x, tumors.30x)){
  to.plot <- c()

  mutation.time.most.likely[[i]] <- c()
  mutation.time.lower[[i]] <- c()
  mutation.time.upper[[i]] <- c()
  
  ## build input matrix
  
  clonal.mutation.matrix.this.tumor <- matrix(unlist(lapply(clonal.mutation.matrix, function(x){x[,i]})), ncol=22,
                                              dimnames = list(rownames(clonal.mutation.matrix$chr1), 
                                                              names(clonal.mutation.matrix)))
  
  segment.length.matrix.this.tumor <- matrix(unlist(lapply(segment.length.matrix, function(x){x[,i]})), ncol=22,
                                             dimnames = list(rownames(clonal.mutation.matrix$chr1), 
                                                             names(clonal.mutation.matrix)))
  
  
  for(j in 1:22){
    
    ## Get mutation counts per copy
    tmp.mut.count <- Mutation.time.converter(clonal.mutation.matrix.this.tumor[,j])
   
    tmp.genome.length <- segment.length.matrix.this.tumor[,j]
    ## Restrict analysis to fragments > 10^7 bp
    tmp.mut.count <- tmp.mut.count[which(tmp.genome.length>10^7)]
    tmp.genome.length <- tmp.genome.length[tmp.genome.length>10^7]

    ## Convert to mutations per haploid genome
    mutation.time.most.likely[[i]] <- c(mutation.time.most.likely[[i]], tmp.mut.count*3.3*10^9/tmp.genome.length)
    mutation.time.lower[[i]] <- c(mutation.time.lower[[i]], 0.5*qchisq(0.025, tmp.mut.count*2)/tmp.genome.length*3.3*10^9)
    mutation.time.upper[[i]] <- c(mutation.time.upper[[i]], 0.5*qchisq(0.975, (tmp.mut.count*2+2))/tmp.genome.length*3.3*10^9)
    
    upscaled.mut.counts <- tmp.mut.count*3.3*10^9/tmp.genome.length
    if(length(upscaled.mut.counts)==0){next}
    names(upscaled.mut.counts) <- paste("chr", j, names(upscaled.mut.counts), sep="_")

    to.plot <- c(to.plot, upscaled.mut.counts)
    names(mutation.time.most.likely[[i]]) <- names(to.plot)
    names(mutation.time.upper[[i]]) <- names(to.plot)
    names(mutation.time.lower[[i]]) <- names(to.plot)
  }
  
  if(length(to.plot)==0){next}

  
  p <- ggplot(data.frame(x=factor(names(mutation.time.most.likely[[i]]), levels = names(mutation.time.most.likely[[i]])[order(mutation.time.most.likely[[i]])]), 
                    y=mutation.time.most.likely[[i]], 
                    ymin=mutation.time.lower[[i]], ymax=mutation.time.upper[[i]]),
         aes(x=x, y=y, ymin=ymin, ymax=ymax)) + geom_col() + geom_errorbar(aes(x=x, ymin=ymin, ymax=ymax)) + 
    scale_y_continuous(limits=c(0, min(max(mutation.time.upper[[i]]), max(mutation.time.upper[[i]]*10))), name = "# Mutations per haploid genome") + 
    theme(axis.text.x = element_text(angle=90), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  print(p)
  
  ## multiply by 2 (diploid genome), scale with 5.1 mutations/day (estimated by Bae et al., 2018) to get time estimates of the time point
  p <- ggplot(data.frame(x=factor(names(mutation.time.most.likely[[i]]), levels = names(mutation.time.most.likely[[i]])[order(mutation.time.most.likely[[i]])]), 
                         y=(mutation.time.most.likely[[i]]*2/5.1+0.1)/7, 
                         ymin=(mutation.time.lower[[i]]*2/5.1+0.1)/7, ymax=(mutation.time.upper[[i]]*2/5.1+0.1)/7),
              aes(x=x, y=y, ymin=ymin, ymax=ymax)) + geom_col() + geom_errorbar(aes(x=x, ymin=ymin, ymax=ymax)) + 
    scale_y_continuous(limits=c(0, min(max(mutation.time.upper[[i]]*2/5.1)/7, max(mutation.time.upper[[i]]*2/5.1/7*10))), name = "Weeks after fertilization") + 
    theme(axis.text.x = element_text(angle=90), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  print(p)
  

}

##############################################################################################################################################
## Estimate the mutational density at ECA and MRCA testing with negative binomial distributions

##############################################################################################################################################
#### Obtain estimates for the mutational density at the tumor's MRCA and test, whether the remaining mutations map to a common ancestor as well. 
## I. For each tumor, determine the mutation time of the lower-order clonal peak (Nmrca). Scale up to entire genome. This should be okay, as the lower-order peak should span most of the genome. Gives us the time point of the MRCA
## II. For each aneuploid fragment, test, whether it conforms to the expectation, which is binomial drawing from Nmrca with f/l, where f is the fragment length and l the total length of the haploid genome
##    Annotate either fragment that doesn't conform
## III. Test for each pre-CNV-fragment whether it conforms to the haploid peak. If so, this fragment is classified as late. If not, it's classified as early
## IV. Compute the total number of mutations among early fragments
## V. Test, whether the remaining of them conform to binomial sampling from the early mutations. 

## do  a confidence interval in the end for the mutation time at mrca and eca, by bootstrapping the segments
mutation.time.mrca <- c()
mutation.time.mrca.lower <- c()
mutation.time.mrca.upper <- c()

earliest.mutation.time <- c()
earliest.mutation.time.lower <- c()
earliest.mutation.time.upper <- c()
gains.at.earliest.time <- list()
gains.not.mapping.to.earliest.time <- list()
other.gains.mapping.to.earliest.time <- list()

mutation.time.eca <- c()
mutation.time.eca.lower <- c()
mutation.time.eca.upper <- c()
gains.at.mrca <- list()
gains.uniquely.mapped.to.eca <- list()
gains.not.maping.to.eca.or.mrca <- list()
monosomic.states.not.matching.mrca <- list()
gains.at.mrca.conforming.eca <- list()

for(i in c(tumors.80x, tumors.30x)){
  print(i)
  
  ## build input matrix
  
  clonal.mutation.matrix.this.tumor <- matrix(unlist(lapply(clonal.mutation.matrix, function(x){x[,i]})), ncol=22,
                                              dimnames = list(rownames(clonal.mutation.matrix$chr1), 
                                                              names(clonal.mutation.matrix)))
  
  segment.length.matrix.this.tumor <- matrix(unlist(lapply(segment.length.matrix, function(x){x[,i]})), ncol=22,
                                              dimnames = list(rownames(clonal.mutation.matrix$chr1), 
                                                              names(clonal.mutation.matrix)))
  
  mrca.eca <- MRCA.ECA.quantification(clonal.mutation.matrix.this.tumor, segment.length.matrix.this.tumor)
  

  mutation.time.mrca[i] <- mrca.eca$mutation.time.mrca
  mutation.time.mrca.lower[i] <- mrca.eca$mutation.time.mrca.lower
  mutation.time.mrca.upper[i] <- mrca.eca$mutation.time.mrca.upper
  
  if(length(mrca.eca$earliest.mutation.time)>0){
    earliest.mutation.time[i] <- mrca.eca$earliest.mutation.time
    earliest.mutation.time.lower[i] <- mrca.eca$earliest.mutation.time.lower
    earliest.mutation.time.upper[i] <- mrca.eca$earliest.mutation.time.upper
  }
  gains.at.earliest.time[[i]] <- mrca.eca$gains.at.earliest.time
  gains.not.mapping.to.earliest.time[[i]] <- mrca.eca$gains.not.mapping.to.earliest.time
  other.gains.mapping.to.earliest.time[[i]] <- mrca.eca$other.gains.mapping.to.earliest.time
  
  mutation.time.eca[i] <- mrca.eca$mutation.time.eca
  mutation.time.eca.lower[i] <- mrca.eca$mutation.time.eca.lower
  mutation.time.eca.upper[i] <- mrca.eca$mutation.time.eca.upper
  gains.at.mrca[[i]] <- mrca.eca$gains.at.mrca
  gains.uniquely.mapped.to.eca[[i]] <- mrca.eca$gains.uniquely.mapped.to.eca
  gains.not.maping.to.eca.or.mrca[[i]] <- mrca.eca$gains.not.maping.to.eca.or.mrca
  monosomic.states.not.matching.mrca[[i]] <- mrca.eca$monosomic.states.not.matching.mrca
  gains.at.mrca.conforming.eca[[i]] <- mrca.eca$gains.at.mrca.conforming.eca
  
}

save(mutation.time.most.likely, mutation.time.upper, mutation.time.lower, earliest.mutation.time, earliest.mutation.time.lower, earliest.mutation.time.upper, gains.at.earliest.time,
     gains.at.mrca, gains.uniquely.mapped.to.eca, gains.not.maping.to.eca.or.mrca, monosomic.states.not.matching.mrca,
     mutation.time.eca, mutation.time.eca.lower, mutation.time.eca.upper, mutation.time.mrca, mutation.time.mrca.lower, mutation.time.mrca.upper,
     gains.at.mrca.conforming.eca, file=paste0(rdata.directory, "MRCA_timing.RData"))

earliest.mutation.time <- earliest.mutation.time[intersect(names(earliest.mutation.time),c(tumors.80x, tumors.30x))]
earliest.mutation.time.lower <- earliest.mutation.time.lower[intersect(names(earliest.mutation.time),c(tumors.80x, tumors.30x))]
earliest.mutation.time.upper <- earliest.mutation.time.upper[intersect(names(earliest.mutation.time),c(tumors.80x, tumors.30x))]
mutation.time.mrca.lower <- mutation.time.mrca.lower[c(tumors.80x, tumors.30x)]
mutation.time.mrca.upper <- mutation.time.mrca.upper[c(tumors.80x, tumors.30x)]
mutation.time.mrca <- mutation.time.mrca[c(tumors.80x, tumors.30x)]
mutation.time.eca <- mutation.time.eca[c(tumors.80x, tumors.30x)]
mutation.time.eca.lower <- mutation.time.eca.lower[c(tumors.80x, tumors.30x)]
mutation.time.eca.upper <- mutation.time.eca.upper[c(tumors.80x, tumors.30x)]
gains.at.earliest.time <- gains.at.earliest.time[intersect(names(gains.at.earliest.time),c(tumors.80x, tumors.30x))]
gains.at.mrca <- gains.at.mrca[intersect(names(gains.at.mrca),c(tumors.80x, tumors.30x))]
gains.uniquely.mapped.to.eca <- gains.uniquely.mapped.to.eca[intersect(names(gains.uniquely.mapped.to.eca),c(tumors.80x, tumors.30x))]
gains.not.maping.to.eca.or.mrca <- gains.not.maping.to.eca.or.mrca[intersect(names(gains.not.maping.to.eca.or.mrca),c(tumors.80x, tumors.30x))]
monosomic.states.not.matching.mrca <- monosomic.states.not.matching.mrca[intersect(names(monosomic.states.not.matching.mrca),c(tumors.80x, tumors.30x))]
gains.at.mrca.conforming.eca <- gains.at.mrca.conforming.eca[intersect(names(gains.at.mrca.conforming.eca),c(tumors.80x, tumors.30x))]

median(mutation.time.mrca)/3.3/10^3
median(mutation.time.eca, na.rm=T)/3.3/10^3
quantile(mutation.time.mrca)/3.3/10^3
quantile(mutation.time.mrca[!is.na(mutation.time.eca)])/3.3/10^3
quantile(mutation.time.eca, na.rm=T)/3.3/10^3

mutation.time.eca[names(earliest.mutation.time)] <- earliest.mutation.time
mutation.time.eca.lower[names(earliest.mutation.time)] <- earliest.mutation.time.lower
mutation.time.eca.upper[names(earliest.mutation.time)] <- earliest.mutation.time.upper
## report this value
quantile(mutation.time.eca, na.rm=T)/3.3/10^3

##############################################################################################################################################
## Alternative check: do segment CIs overlap with MRCA?
evidence.for.eca.from.individual.ci <- c()

for(i in 1:length(c(tumors.80x, tumors.30x))){
  fragments.timing.mrca <- mutation.time.upper[[c(tumors.80x, tumors.30x)[i]]][sapply(names(mutation.time.upper[[c(tumors.80x, tumors.30x)[i]]]), function(x){
    x <- strsplit(x, split="_")[[1]][4]
    if(is.na(x) || x=="I"){T}else{
      F
    }
  })]
  evidence.for.eca.from.individual.ci[i] <- any(mutation.time.upper[[i]] < min(fragments.timing.mrca))
}
names(evidence.for.eca.from.individual.ci) <- c(tumors.80x, tumors.30x)

intersect(names(evidence.for.eca.from.individual.ci==T), names(mutation.time.eca[!is.na(mutation.time.eca)]))
setdiff(names(evidence.for.eca.from.individual.ci==T), names(mutation.time.eca[!is.na(mutation.time.eca)]))
setdiff(names(mutation.time.eca[!is.na(mutation.time.eca)]), names(evidence.for.eca.from.individual.ci==T))
## would increase number of tumors with ECA

##############################################################################################################################################
### Extract for each tumor the VAF distribution for each ploidy state. Exclude sex chromosomes

vafs.all.tumors <- list(0)
genome.size.all.tumors <- list(0)

load(paste0(rdata.directory, "Clonal_mutations_different_ploidies.RData"))

purities.all.tumors <- purity
ploidies.all.tumors <- ploidy
for(i in c(tumors.80x, tumors.30x)){
  
  print(i)
  
  if(i %in% tumors.80x){
    data.directory <- data.directory.80x
  }else{
    data.directory <- data.directory.30x
  }
  
  ## Read in copy number information and extract ploidy/purity, as before
  aceseq <- list.files(paste0(data.directory, "/", i, "/", cnv.directory), pattern="comb_pro_extra")[1]
  while(is.na(aceseq)){
    aceseq <- list.files(paste0(data.directory, "/", i, "/", cnv.directory), pattern="comb_pro_extra")[1]
    print(i)
  }

  ## read in the mutation file
  files <- list.files(paste0(data.directory, "/", i, "/", snv.directory, "/"), pattern="somatic_snvs_conf_8_to_10")[1]
  mutations <- read.vcf(paste0(data.directory, "/", i, "/", snv.directory, "/", files))
  mutations$vcf <- mutations$vcf[!mutations$vcf$CHROM %in% c("X", "Y"), ]
  
  purity <- purities.all.tumors[i]
  ploidy <- ploidies.all.tumors[i]
  
  copy.number.info <- read.delim(file=paste0(data.directory, "/", i, "/", cnv.directory, "/", aceseq), sep="\t", stringsAsFactors = F)
  copy.number.info <- copy.number.info[!copy.number.info$X.chromosome %in% c("X", "Y"),]
  ## obtain the coverage ratios for the mutations of interest
  
  ## Extract copy number info for each mutation
  cnv.info.per.mutation <- Extract.copy.number.info.per.SSNV(mutations, copy.number.info)
  ## obtain the coverage ratios at mutated loci
  coverage.ratios <- cnv.info.per.mutation$coverage.ratio
  bafs <- cnv.info.per.mutation$baf
  genotype <- cnv.info.per.mutation$genotype
  tcn <- cnv.info.per.mutation$tcn


  ## Extract readcounts of reference and alternative bases
  readcounts <- Extract.info.from.vcf(mutations, info="readcounts")
  
  
    #######################################################################
    ## Iterate through all copy number states
  
  vafs.this.tumor <- list(0)
  genome.size.this.tumor <- list(0)
    

  ## Plot separately for each copy number
    for(k in unique(copy.number.indicator)){
      
      expected.coverage.ratio <- (k*purity + (1-purity)*2)/(ploidy*purity+(1-purity)*2)
      
      readcounts. <- readcounts[((coverage.ratios>(expected.coverage.ratio-0.1) & coverage.ratios<(expected.coverage.ratio+0.1) & !is.na(coverage.ratios)) |
                                   (tcn ==k & !is.na(tcn))) ,,drop=F]
      
      
      vafs.this.tumor[[k]] <- readcounts.
      
      
      genome.size <- sum(as.numeric(copy.number.info[(copy.number.info$tcnMeanRaw>(expected.coverage.ratio-0.1) & copy.number.info$tcnMeanRaw<(expected.coverage.ratio+0.1)) |
                                                       (as.numeric(copy.number.info$TCN)==k & !is.na(as.numeric(copy.number.info$TCN))) ,]$end)-
                           as.numeric(copy.number.info[(copy.number.info$tcnMeanRaw>(expected.coverage.ratio-0.1) & copy.number.info$tcnMeanRaw<(expected.coverage.ratio+0.1)) |
                                                         (as.numeric(copy.number.info$TCN)==k & !is.na(as.numeric(copy.number.info$TCN))),]$start))
      
      genome.size.this.tumor[[k]] <- genome.size
      
      
      for(l in B.allele.indicator[which(copy.number.indicator==k)]){
      
      
        readcounts. <- readcounts[((coverage.ratios>(expected.coverage.ratio-0.1) & coverage.ratios<(expected.coverage.ratio+0.1) & !is.na(coverage.ratios)) |
                                     (tcn ==k & !is.na(tcn))) &
                                    (((bafs < (max(l/k, 1-l/k)+0.05) & bafs > (max(l/k, 1-l/k)-0.05)) & !is.na(bafs) | (is.na(bafs) & l==k/2)) |
                                       (genotype==paste(k-l, l, sep=":") & !is.na(genotype)) | 
                                       (genotype==paste(l, k-l, sep=":") & !is.na(genotype))),,drop=F]
        
        if(nrow(readcounts.)==0){next}
      
        

        prob.clonal <- l*purity/(purity*k + (1-purity)*2)
        monosomic.prob.clonal <- purity/(purity*k + (1-purity)*2)
        

       p <-  ggplot(data.frame(VAF=readcounts.[,2]/rowSums(readcounts.)), aes(x=VAF)) + geom_histogram(binwidth = 0.01) + 
          geom_vline(xintercept = prob.clonal, col="firebrick", linetype=2) + geom_vline(xintercept=monosomic.prob.clonal, col="firebrick", linetype=2)+ 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(name="# Mutations") +
          scale_x_continuous(limits=c(0,1)) + ggtitle(paste0("Copy number = ", k, "B allele = ", l))
       
       print(p)

      
      }
    }

  
  vafs.all.tumors[[i]] <- vafs.this.tumor
  genome.size.all.tumors[[i]] <- genome.size.this.tumor
  save(vafs.all.tumors, genome.size.all.tumors, file=paste0(rdata.directory, "Vafs_all_tumors.RData"))
  
}

save(vafs.all.tumors, genome.size.all.tumors, file=paste0(rdata.directory, "Vafs_all_tumors.RData"))



