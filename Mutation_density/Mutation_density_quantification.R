##############################################################################################################################################
## Quantify the densitiy of amplified and non-amplified clonal mutations
##############################################################################################################################################
## load libraries
library(Hmisc)
library(scales)
library(mixtools)
library(bedr)

## Required input: SNVs, Indels, ACEseq output (*comb_pro_extra*). The input files should be organized per tumor
## For each tumor store SNVs and the ACEseq output in subfolders named SNVs, ACEseq

## fill in main directory
data.directory <- ""
## specify tumors
tumors <- c("Example_tumor1_primary", "Example_tumor1_relapse", "Example_tumor2_primary",
            "Example_tumor3_primary", "Example_tumor3_relapse")
##  primary tumor of tumor pairs
primary.of.tumor.pairs <- c("Example_tumor1_primary", "Example_tumor3_primary")
## and relapse tumor of tumor pairs
relapse.of.tumor.pairs <- c("Example_tumor1_relapse", "Example_tumor3_relapse")

source("Estimate_sampling_bias_from_tumor_pairs.R")

##############################################################################################################################################
## For each tumor, we want to know the number of clonal mutations that lie on one or multiple copies of a gained chromosome. 
## In order to compute densities, we also need to report the fraction of the genome at a specific copy number.

## Store ACEseq estimates of purity and ploidy for each tumor
purity <- rep(0, length(tumors))
ploidy <- purity

## Store the clonal mutations present on any number of alleles and store this information for each chromosome separately in a list.
## The roman numbers correspond to the number of alleles carrying a mutation.

clonal.mutation.matrix. <- matrix(0, nrow=1 + 2 + 3 + 4 , ncol=length(tumors),
                                 dimnames=list(c("monosomic", "disomic_I", "disomic_II", "trisomic_I", "trisomic_II", "trisomic_III",
                                                 "tetrasomic_I", "tetrasomic_II", "tetrasomic_III", "tetrasomic_IV"),
                                               tumors))
## one matrix per autosome
clonal.mutation.matrix <- list()

for(i in 1:22){
  clonal.mutation.matrix[[i]] <- clonal.mutation.matrix.
  names(clonal.mutation.matrix)[i] <- paste0("chr",i)
}
rm(clonal.mutation.matrix.)


## analogously, we generate a table to store the genomic length associated with each copy number state

genome.length.matrix <- clonal.mutation.matrix

## Indicator vectors for the copy number and the B allele count; all combinations for CN <= 4

copy.number.indicator <- c(unlist(sapply(1:4, function(x){rep(x, each=x)})))
B.allele.indicator <- c(unlist(sapply(1:4, function(x){1:x})))


## Iterate through the tumors, iterate through the autosomes and iterate through the different B allele frequencies

## First loop: tumors
for(i in tumors){
  
  print(i)
  
  ## Read in purity-ploidy information (from Aceseq and validated with FACS information where available):
  aceseq <- list.files(paste0(data.directory, "/", i, "/ACEseq/"), pattern="comb_pro_extra")[1]
  if(is.na(aceseq)){next}
  copy.number.info. <- read.delim(file=paste0(data.directory, "/", i, "/ACEseq/", aceseq), sep="\t", stringsAsFactors = F)
  
  ## read in the mutation file
  files <- list.files(paste0(data.directory, "/", i, "/SNVs/"), pattern="somatic_snvs_conf_8_to_10")[1]
  
  mutations <- read.vcf(paste0(data.directory, "/", i, "/SNVs/", files))
  
  ## obtain the coverage ratios at mutated loci
  coverage.ratios. <- apply(mutations$vcf, 1, function(x){
    x <- unlist(x)
    tmp <- copy.number.info.[copy.number.info.$X.chromosome==x[1],]
    tmp <- tmp[as.numeric(tmp$start) <= as.numeric(x[2]) & as.numeric(tmp$end) >= as.numeric(x[2]),]
    if(nrow(tmp)>0){
      unlist(tmp[28])
    }else{
      NA
    }
  })
  ## and the B-allele frequencies
  bafs. <- apply(mutations$vcf, 1, function(x){
    x <- unlist(x)
    tmp <- copy.number.info.[copy.number.info.$X.chromosome==x[1],]
    tmp <- tmp[as.numeric(tmp$start) <= as.numeric(x[2]) & as.numeric(tmp$end) >= as.numeric(x[2]),]
    if(nrow(tmp)>0){
      unlist(tmp[30])
    }else{
      NA
    }
  })
  ## and the genotype
  genotype. <- apply(mutations$vcf, 1, function(x){
    x <- unlist(x)
    tmp <- copy.number.info.[copy.number.info.$X.chromosome==x[1],]
    tmp <- tmp[as.numeric(tmp$start) <= as.numeric(x[2]) & as.numeric(tmp$end) >= as.numeric(x[2]),]
    if(nrow(tmp)>0){
      unlist(tmp[36])
    }else{
      NA
    }
  })
  ## and the suggested copy number
  tcn. <- apply(mutations$vcf, 1, function(x){
    x <- unlist(x)
    tmp <- copy.number.info.[copy.number.info.$X.chromosome==x[1],]
    tmp <- tmp[as.numeric(tmp$start) <= as.numeric(x[2]) & as.numeric(tmp$end) >= as.numeric(x[2]),]
    if(nrow(tmp)>0){
      as.numeric(unlist(tmp[37]))
    }else{
      NA
    }
  })
  
  ## Extract ploidy and purity estimates from ACEseq output
  ploidy. <- sapply(aceseq, function(x){strsplit(x, split="extra")[[1]][2]})
  purity. <- sapply(ploidy., function(x){strsplit(x, split="_")[[1]][2]})
  purity. <- as.numeric(sapply(purity., function(x){strsplit(x, split=".txt")[[1]][1]}))
  ploidy. <- as.numeric(sapply(ploidy., function(x){strsplit(x, split="_")[[1]][1]}))

  purity[i] <- purity.
  ploidy[i] <- ploidy.
  

  ## Extract readcounts of reference and alternative bases
  readcounts. <- t(sapply(mutations$vcf$INFO, function(x){
    x <- strsplit(x, split=";")[[1]][6]
    x <- strsplit(x, split="=")[[1]][2]
    x <- as.numeric(strsplit(x, split=",")[[1]])
    x <- unname(x)
    c(sum(x[c(1,2)]), sum(x[c(3,4)]))
  }))
  rownames(readcounts.) <- paste(mutations$vcf$GENE, mutations$vcf$POS, sep=".")
  colnames(readcounts.) <- c("REF", "ALT")
  
  ## Correct purity estimate for tumors NBE26, NBE40 manually
  if(i %in% c("NBE40", "NBE26")){
    mixmdl = normalmixEM(readcounts.[!is.na(rowSums(readcounts.)),2]/rowSums(readcounts.[!is.na(rowSums(readcounts.)),]))
    purity. <- max(mixmdl$mu)*2
    while(purity.>1){
      mixmdl = normalmixEM(readcounts.[!is.na(rowSums(readcounts.)),2]/rowSums(readcounts.[!is.na(rowSums(readcounts.)),]))
      purity. <- max(mixmdl$mu)*2
    }

    purity[i] <- purity.
  }
  

  #######################################################################
  ## Iterate through all autosomes
  
  for(j in 1:22){
    
    ## Extract copy number info etc. for the current chromosome
    copy.number.info <- copy.number.info.[copy.number.info.$X.chromosome==j,,drop=F]
    if(nrow(copy.number.info)==0){next}
    bafs <- bafs.[which(mutations$vcf$CHROM==j)]
    if(length(bafs)==0){next}
    coverage.ratios <- coverage.ratios.[which(mutations$vcf$CHROM==j)]
    tcn <- tcn.[which(mutations$vcf$CHROM==j)]
    genotype <- genotype.[which(mutations$vcf$CHROM==j)]
    if(length(coverage.ratios)==0){next}
    readcounts <- readcounts.[which(mutations$vcf$CHROM==j),,drop=F]
    if(nrow(readcounts)==0){next}
    
    #######################################################################
    ## Iterate through all copy number states
    
    for(k in unique(copy.number.indicator)){
      
      
      #######################################################################
      ## Iterate through all A/B combinations. The purpose here is now to determine the number of mutations in the clonal peaks.
      ## In regions that are monosomic or disomic there is only one clonal peak, which does not carry much information.
      ## However, in higher orders, we can use the number of mutations on more than one allele as a surrogate for the time that has passed since conception
      ## and the amplification event. 
      ## We fit all clonal peaks (for the lower-order peak only the right hand side to avoid contamination with subclonal mutations) to the sum of two binomials in dependence of the number of mutations belonging to 
      ## each peak. 
      
      for(l in B.allele.indicator[which(copy.number.indicator==k)][ceiling(length(B.allele.indicator[which(copy.number.indicator==k)])/2):
                                                                   length(B.allele.indicator[which(copy.number.indicator==k)])]){
        

        ## Frequency of the clonal peak associated with mutations on monosomic chromosomes
        monosomic.prob.clonal <- purity[i]/(purity[i]*k + (1-purity[i])*2)
        
        ## Expected frequencies of lower- and higher-order clonal peaks
        clonal.peaks <- sort(unique(sapply(c(1, l, k-l), function(x){x*purity[i]/(purity[i]*k + (1-purity[i])*2)})))
        clonal.peaks <- clonal.peaks[clonal.peaks!=0]
        
        ## Assign regions to a specific copy number with tolerance +/- 0.1 from expected coverage ratio
        expected.coverage.ratio <- (k*purity[i] + (1-purity[i])*2)/(ploidy[i]*purity[i]+(1-purity[i])*2)
        genome.length.matrix[[j]][which(copy.number.indicator==k & B.allele.indicator %in% c(k-l,l)),i] <-  genome.length.matrix[[j]][which(copy.number.indicator==k & B.allele.indicator %in% c(k-l,l)),i] +
          rep(sum(as.numeric(copy.number.info[((copy.number.info$tcnMeanRaw>(expected.coverage.ratio-0.1) & copy.number.info$tcnMeanRaw<(expected.coverage.ratio+0.1))|
                                                                                                                                   (as.numeric(copy.number.info$TCN)==k & !is.na(as.numeric(copy.number.info$TCN))) ) &
                                                                                                                                  (((copy.number.info$BAF < (max(l/k, 1-l/k)+0.05) & copy.number.info$BAF > (max(l/k, 1-l/k)-0.05)) &
                                                                                                                                      !is.na(copy.number.info$BAF) | (is.na(copy.number.info$BAF) & l==k/2)) |
                                                                                                                                     ((copy.number.info$genotype==paste(k-l, l, sep=":") | copy.number.info$genotype==paste(l, k-l, sep=":")) &
                                                                                                                                        !is.na(copy.number.info$genotype))),]$end)-
                                                                                                      as.numeric(copy.number.info[((copy.number.info$tcnMeanRaw>(expected.coverage.ratio-0.1) & copy.number.info$tcnMeanRaw<(expected.coverage.ratio+0.1))|
                                                                                                                                     (as.numeric(copy.number.info$TCN)==k & !is.na(as.numeric(copy.number.info$TCN))) ) &
                                                                                                                                    (((copy.number.info$BAF < (max(l/k, 1-l/k)+0.05) & copy.number.info$BAF > (max(l/k, 1-l/k)-0.05)) &
                                                                                                                                        !is.na(copy.number.info$BAF) | (is.na(copy.number.info$BAF) & l==k/2)) |
                                                                                                                                       ((copy.number.info$genotype==paste(k-l, l, sep=":") | copy.number.info$genotype==paste(l, k-l, sep=":")) &
                                                                                                                                          !is.na(copy.number.info$genotype))),]$start)), ifelse(k==l, 1, length(unique(c(l, k-l)))))
        
        
        
      
        ## the monosomic fraction is putatively measured for several states, thus add up
        if(l!=1 & (k-l)!=1){
          genome.length.matrix[[j]][which(copy.number.indicator==k & B.allele.indicator==1),i] <- genome.length.matrix[[j]][which(copy.number.indicator==k & B.allele.indicator==1),i] +
            sum(as.numeric(copy.number.info[((copy.number.info$tcnMeanRaw>(expected.coverage.ratio-0.1) & copy.number.info$tcnMeanRaw<(expected.coverage.ratio+0.1))|
                                               (as.numeric(copy.number.info$TCN)==k & !is.na(as.numeric(copy.number.info$TCN))) ) &
                                              (((copy.number.info$BAF < (max(l/k, 1-l/k)+0.05) & copy.number.info$BAF > (max(l/k, 1-l/k)-0.05)) &
                                                  !is.na(copy.number.info$BAF) | (is.na(copy.number.info$BAF) & l==k/2)) |
                                                 ((copy.number.info$genotype==paste(k-l, l, sep=":") | copy.number.info$genotype==paste(l, k-l, sep=":")) &
                                                    !is.na(copy.number.info$genotype))),]$end)-
                  as.numeric(copy.number.info[((copy.number.info$tcnMeanRaw>(expected.coverage.ratio-0.1) & copy.number.info$tcnMeanRaw<(expected.coverage.ratio+0.1))|
                                                 (as.numeric(copy.number.info$TCN)==k & !is.na(as.numeric(copy.number.info$TCN))) ) &
                                                (((copy.number.info$BAF < (max(l/k, 1-l/k)+0.05) & copy.number.info$BAF > (max(l/k, 1-l/k)-0.05)) &
                                                    !is.na(copy.number.info$BAF) | (is.na(copy.number.info$BAF) & l==k/2)) |
                                                   ((copy.number.info$genotype==paste(k-l, l, sep=":") | copy.number.info$genotype==paste(l, k-l, sep=":")) &
                                                      !is.na(copy.number.info$genotype))),]$start))
          
          
        }
     
        ## If the current copy number state does not exist on the present autosome, continue
        if(genome.length.matrix[[j]][which(copy.number.indicator==k & B.allele.indicator==l),i]==0){next}
        
        ## Count the number of mutations that lie on the respective copy numbers
        ## 1st: choose all mutations that lie on loci with current copy number
        mutations.at.copy.number <-  readcounts[((coverage.ratios>(expected.coverage.ratio-0.1) & coverage.ratios < (expected.coverage.ratio+0.1) & !is.na(coverage.ratios)) |
                                                                             tcn==k & !is.na(tcn)) &
                                                                            (((bafs < (max(l/k, 1-l/k)+0.05) & bafs > (max(l/k, 1-l/k)-0.05)) & !is.na(bafs) | (is.na(bafs) & l==k/2)) |
                                                                               genotype==paste(k-l, l, sep=":") | genotype==paste(l, k-l, sep=":")),,drop=F]
        
        if(nrow(mutations.at.copy.number)==0){next}
        

        ## 2nd we now count all mutations that lie on the right hand side of the monosomic clonal peak to remove subclonal mutations

        mutations.at.copy.number <- mutations.at.copy.number[mutations.at.copy.number[,2]/rowSums(mutations.at.copy.number)>=monosomic.prob.clonal,,drop=F]
        mutations.at.copy.number <- mutations.at.copy.number[!is.na(rowSums(mutations.at.copy.number)),,drop=F]
        
        if(nrow(mutations.at.copy.number)==0){next}
          
        ## If there's only one mutation, assign it to its most likely state by comparing its VAF with the clonal frequencies
        if(nrow(mutations.at.copy.number) == 1){
          which.peak <- which.min(c(clonal.peaks - mutations.at.copy.number[,2]/rowSums(mutations.at.copy.number))^2)
          which.peak <- c(1, sort(c(l, k-l)))[which.peak]
          
          clonal.mutation.matrix[[j]][which(copy.number.indicator==k & B.allele.indicator==which.peak),i] <- 1
          next
        }
        
        ## In case of monosomic or heterozygous disomic regions, there is only one clonal peak. Thus  assign all mutations to that peak
        if(k %in% c(1,2) & l==1){
          clonal.mutation.matrix[[j]][which(copy.number.indicator==k & B.allele.indicator==l),i] <- nrow(mutations.at.copy.number) 
        next}
        
        
        ## For the remaining cases, estimate the sizes of the clonal peaks using a binomial mixture model.
        ## Define the posterior probability of a mixing model with mixing factor p. Scan p over 0 - 1 and select the most likely one.
        
        p.priors <- seq(0, 1, 0.01)
        posteriors <- sapply(p.priors, function(p){
          sum(apply(mutations.at.copy.number, 1, function(x){
            L <-  dbinom(x = x[2], size=sum(x), prob = clonal.peaks)
            P <- L/sum(L)
            log(sum(c(p, rep(1-p, length(clonal.peaks) - 1))*P))
          }))
        })
        
        
        clonal.mutation.matrix[[j]][which(copy.number.indicator==k & B.allele.indicator %in% c(1, l, k-l)),i] <- c(nrow(mutations.at.copy.number)*p.priors[which.max(posteriors)],
                                                                                                                   rep(nrow(mutations.at.copy.number)*(1-p.priors[which.max(posteriors)]),
                                                                                                                       length(clonal.peaks)-1)) +
          clonal.mutation.matrix[[j]][which(copy.number.indicator==k & B.allele.indicator %in% c(1, l, k-l)),i] 
      
              
      }
    }
  }
}

##############################################################################################################################################
## Compute the mutational density at each genomic fragment 

## To this end, the monosomic peak needs to be multiplied with two, since we only estimated the RHS. 

# and store the mutation time of the 1st event

mutation.time.most.likely <- list()

for(i in tumors){
  mutation.time.most.likely[[i]] <- c()

  for(j in 1:22){
    
    ## The lower-order clonal peak was only quantified on its RHS. Thus multiply by 2. Moreover, if both alleles have the same count, divide by 2, 
    ## accounting for the 2 alleles (e.g. 2:2, 3:3)
    tmp.mut.count <- clonal.mutation.matrix[[j]][,i]*c(2, 2, 1, 2, 1, 1, 2, 1/2, 1, 1)
    ## In addition, add the number of mutations, which were lost from the lower-order peak due to amplification, to each copy of an amplified allele
    tmp.mut.count[c(1,2,4,7)] <- tmp.mut.count[c(1,2,4,7)] + c(0,
                                                                           sum(tmp.mut.count[c(3)]*2), 
                                                                           sum(tmp.mut.count[c(5,6)]*c(2,3)) , 
                                                                           sum(tmp.mut.count[c(8,9,10)]*c(4,3,4)))
    
    ## divide by the copy number to obtain mutations per copy
    tmp.mut.count[c(1,2,4,7)] <- tmp.mut.count[c(1,2,4,7)]/(c(1,2,3,4))
    
    tmp.genome.length <- genome.length.matrix[[j]][,i]
    ## Restrict analysis to fragments > 10^7 bp
    tmp.mut.count <- tmp.mut.count[which(tmp.genome.length>10^7)]
    tmp.genome.length <- tmp.genome.length[tmp.genome.length>10^7]

    if(length(tmp.mut.count)==0){next}
    ## Convert to mutations per haploid genome
    names(tmp.mut.count) <- paste("chr", j, names(tmp.mut.count), sep="_")
    mutation.time.most.likely[[i]] <- c(mutation.time.most.likely[[i]], tmp.mut.count*3.3*10^9/tmp.genome.length)

  }

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

## Estimated mutation time at mrca (store for each tumor)
mutation.time.mrca <- c()
mutation.time.mrca.lower <- c()
mutation.time.mrca.upper <- c()

## if there are multiple detectable time points prior to the MRCA, store the earliest time point here
earliest.mutation.time <- c()
earliest.mutation.time.lower <- c()
earliest.mutation.time.upper <- c()

## Estimated mutation time at ECA; store for each tumor
mutation.time.eca <- c()
mutation.time.eca.lower <- c()
mutation.time.eca.upper <- c()

## Store aneuploidies that can be timed at the MRCA (each list entry is one tumor)
aneuploidies.at.mrca <- list()
## Store  aneuploidies that cannot be timed at the MRCA but match the ECA
remaining.aneuploidies.conforming.single.event <- list()
## Store  aneuploidies that cannot be timed at the MRCA nor at the ECA
remaining.aneuploidies.not.conforming.single.event <- list()
haploids.not.matching.mrca <- list()
## Store aneuploidies that match both, MRCA and ECA
aneuploidies.at.mrca.conforming.eca <- list()

## Aneuplodies acquired at the earliest time point
aneuploidies.at.earliest.time <- list()
remaining.aneuploidies.not.conforming.earliest.event <- list()
other.aneuploidies.conforming.earliest.time <- list()

for(i in tumors){
  print(i)

  mrca.mutation.count <- c()
  mrca.genome.length <- c()
  
  for(j in 1:22){ 
    ## The lower-order clonal peak was quantified on its RHS only. Thus multiply by 2. Moreover, if both alleles have the same count, divide by 2, 
    ## accounting for the 2 alleles (e.g. 2:2, 3:3)
    tmp.mut.count <- clonal.mutation.matrix[[j]][,i]*c(2, 2, 1, 2, 1, 1, 2, 1/2, 1, 1)
    
    ## In addition, add the number of mutations, which do not contriubte to the monosomic peak due to amplification, on the respective number of alleles
    tmp.mut.count[c(1,2,4,7)] <- tmp.mut.count[c(1,2,4,7)] + c(0,
                                                               sum(tmp.mut.count[c(3)]*2), 
                                                               sum(tmp.mut.count[c(5,6)]*c(2,3)) , 
                                                               ## times 4 since we already accounted for the 2 alleles above
                                                               sum(tmp.mut.count[c(8,9,10)]*c(4,3,4)) )
    ## divide by the copy number to obtain mutations per copy
    tmp.mut.count[c(1,2,4,7)] <- tmp.mut.count[c(1,2,4,7)]/(c(1,2,3,4))
    
    tmp.genome.length <- genome.length.matrix[[j]][,i]
    ## Take only lower-order peaks
    tmp.mut.count <- tmp.mut.count[c(1,2,4,7)]
    tmp.genome.length <- tmp.genome.length[c(1,2,4,7)]
    
    ## Restrict analysis to fragments > 10^7 bp
    tmp.mut.count <- tmp.mut.count[which(tmp.genome.length>10^7)]
    tmp.genome.length <- tmp.genome.length[tmp.genome.length>10^7]
    
    if(length(tmp.genome.length)==0){next}
    
    names(tmp.genome.length) <- paste("chr", j, names(tmp.genome.length), sep="_")
    
    mrca.mutation.count <- c(mrca.mutation.count, tmp.mut.count)
    mrca.genome.length <- c(mrca.genome.length, tmp.genome.length)
    
  }
  if(length(mrca.mutation.count)==0){next}
  
  names(mrca.mutation.count) <- names(mrca.genome.length)
  
  ## subtract the estimated error of clonal mutations due to a sampling bias
  mutation.time.mrca <- c(mutation.time.mrca, sum(mrca.mutation.count)*3.3*10^9/sum(mrca.genome.length)*(1- mean(clonal.mutations.false.positives/clonal.mutations.all)))
  names(mutation.time.mrca)[length(mutation.time.mrca)] <- i
  
  ## bootstrap upper and lower limits of the mutation time; + subtract the number of false positive mutations 
  bootstrapped.time <- sapply(1:1000, function(x){
    res <- sample(x=1:length(mrca.mutation.count), size=length(mrca.mutation.count), prob= mrca.genome.length, replace=T)
    res <- sum(mrca.mutation.count[res])*3.3*10^9/sum(mrca.genome.length[res])
    res <- res - res*rnorm(n=1, mean=mean(clonal.mutations.false.positives/clonal.mutations.all), sd=sd(clonal.mutations.false.positives/clonal.mutations.all))
  })
  mutation.time.mrca.lower <- c(mutation.time.mrca.lower,quantile(bootstrapped.time, 0.025))
  mutation.time.mrca.upper <- c(mutation.time.mrca.upper,quantile(bootstrapped.time, 0.975))
  names(mutation.time.mrca.lower)[length(mutation.time.mrca.lower)] <- i
  names(mutation.time.mrca.upper)[length(mutation.time.mrca.upper)] <- i
  
  ## now, test whether the individual fragments match to a joint time point. Test with a negative binomial distribution to account for overdispersion of the
  ## local mutation rate: 
  
  does.fragment.match <- apply(rbind(mrca.mutation.count, mrca.genome.length), 2, function(x){
    if(x[1] < round(sum(mrca.mutation.count))* x[2]/sum(mrca.genome.length)){
      test <- pnbinom(q = x[1], size = round(sum(mrca.mutation.count)), prob = round(sum(mrca.mutation.count))/(round(sum(mrca.mutation.count))*(1+x[2]/sum(mrca.genome.length))))
    }else{
      test <- pnbinom(q = x[1], size = round(sum(mrca.mutation.count)), prob = round(sum(mrca.mutation.count))/(round(sum(mrca.mutation.count))*(1+x[2]/sum(mrca.genome.length))), lower.tail = F)
    }
    test})
  haploid.fragment.names <- names(does.fragment.match)
  ## collect the p-values and adjust them later, once contributions from early clonal mutations are accounted for, too
  p.values <- does.fragment.match
  
  ## now, test whether mutation counts on amplified fragments are consistent with the mrca
  mean.fp <- mean(clonal.mutations.false.positives/clonal.mutations.all)
  eca.mutation.counts <- c()
  eca.genome.length <- c()
  eca.mutation.counts.at.mrca <- c()
  eca.genome.length.at.mrca <- c()
  aneuploid.fragment.name <- c()
  aneuploid.mutation.counts <- c()
  aneuploid.genome.lengths <- c()
  for(j in 1:22){
    tmp.mut.count <- clonal.mutation.matrix[[j]][,i]*c(2, 2, 1, 2, 1, 1, 2, 1/2, 1, 1)
    tmp.genome.length <- genome.length.matrix[[j]][,i]
    
    tmp.mut.count <- tmp.mut.count[-c(1,2,4,7)]
    tmp.genome.length <- tmp.genome.length[-c(1,2,4,7)]
    
    ## Restrict analysis to fragments > 10^7 bp
    tmp.mut.count <- tmp.mut.count[which(tmp.genome.length>10^7)]
    tmp.genome.length <- tmp.genome.length[tmp.genome.length>10^7]
    
    if(length(tmp.genome.length)==0){next}
    
    names(tmp.genome.length) <- paste("chr", j, names(tmp.genome.length), sep="_")
    names(tmp.mut.count) <- names(tmp.genome.length)
    
    ## does the fragment conform to the mrca? Test with negative binomial distribution
    
    does.fragment.match <- apply(rbind(tmp.mut.count, tmp.genome.length), 2, function(x){
      if(x[1] <= round(sum(mrca.mutation.count))*(1-mean.fp)* x[2]/sum(mrca.genome.length)){
        test <- pnbinom(q = x[1], size = round(sum(mrca.mutation.count)*(1-mean.fp)), prob = round(sum(mrca.mutation.count)*(1-mean.fp))/(round(sum(mrca.mutation.count)*(1-mean.fp))*(1+x[2]/sum(mrca.genome.length))))
      }else{
        test <- pnbinom(q = x[1], size = round(sum(mrca.mutation.count)*(1-mean.fp)), prob = round(sum(mrca.mutation.count)*(1-mean.fp))/(round(sum(mrca.mutation.count)*(1-mean.fp))*(1+x[2]/sum(mrca.genome.length))), lower.tail = F)
      }
      test})
    aneuploid.fragment.name <- c(aneuploid.fragment.name, names(does.fragment.match))
    p.values <- c(p.values, does.fragment.match)
    
    aneuploid.genome.lengths <- c(aneuploid.genome.lengths, tmp.genome.length)
    aneuploid.mutation.counts <- c(aneuploid.mutation.counts, tmp.mut.count)
    
  }
  
  ## Adjust p values and cut off at 0.01. Fragments with p < 0.01 likely arose from a different time point
  adjusted.p.values <- p.adjust(p.values)
  haploids.not.matching.mrca[[i]] <- names(adjusted.p.values[haploid.fragment.names][adjusted.p.values[haploid.fragment.names]<0.01])
  aneuploidies.at.mrca[[i]] <- names(adjusted.p.values[aneuploid.fragment.name][adjusted.p.values[aneuploid.fragment.name]>=0.01])
  eca.mutation.counts.at.mrca <- aneuploid.mutation.counts[names(adjusted.p.values[aneuploid.fragment.name][adjusted.p.values[aneuploid.fragment.name]>=0.01])]
  eca.genome.length.at.mrca <- aneuploid.genome.lengths[names(adjusted.p.values[aneuploid.fragment.name][adjusted.p.values[aneuploid.fragment.name]>=0.01])]
  eca.mutation.counts <- aneuploid.mutation.counts[names(adjusted.p.values[aneuploid.fragment.name][adjusted.p.values[aneuploid.fragment.name]<0.01])]
  eca.genome.length <- aneuploid.genome.lengths[names(adjusted.p.values[aneuploid.fragment.name][adjusted.p.values[aneuploid.fragment.name]<0.01])]
  
  ## From the non-mapping fragments estimate a time of origin and test whether all of them conform to it.
  ## here, fragments with counts higher than the mrca should be excluded. OK as these are likely outliers
  not.used.for.quantification <- which(eca.mutation.counts*3.3*10^9/eca.genome.length > mutation.time.mrca[length(mutation.time.mrca)])
  if(length(not.used.for.quantification) > 0){
    mutation.time.eca <- c(mutation.time.eca, sum(eca.mutation.counts[-not.used.for.quantification])*3.3*10^9/sum(eca.genome.length[-not.used.for.quantification]))
    ## bootstrap upper and lower limits of the mutation time
    if(length(eca.mutation.counts[-not.used.for.quantification])>0){
      bootstrapped.time <- sapply(1:1000, function(x){
        res <- sample(x=1:length(eca.mutation.counts[-not.used.for.quantification]), size=length(eca.mutation.counts[-not.used.for.quantification]), prob= eca.genome.length[-not.used.for.quantification], replace=T)
        sum(eca.mutation.counts[-not.used.for.quantification][res])*3.3*10^9/sum(eca.genome.length[-not.used.for.quantification][res])
      })
      mutation.time.eca.lower <- c(mutation.time.eca.lower,quantile(bootstrapped.time, 0.025))
      mutation.time.eca.upper <- c(mutation.time.eca.upper,quantile(bootstrapped.time, 0.975))
      names(mutation.time.eca.lower)[length(mutation.time.eca.lower)] <- i
      names(mutation.time.eca.upper)[length(mutation.time.eca.upper)] <- i
    }else{
      mutation.time.eca.lower <- c(mutation.time.eca.lower,NA)
      mutation.time.eca.upper <- c(mutation.time.eca.upper,NA)
      names(mutation.time.eca.lower)[length(mutation.time.eca.lower)] <- i
      names(mutation.time.eca.upper)[length(mutation.time.eca.upper)] <- i
    }
  }else{
    mutation.time.eca <- c(mutation.time.eca, sum(eca.mutation.counts)*3.3*10^9/sum(eca.genome.length))
    ## bootstrap upper and lower limits of the mutation time
    if(length(eca.mutation.counts)>0){
      bootstrapped.time <- sapply(1:1000, function(x){
        res <- sample(x=1:length(eca.mutation.counts), size=length(eca.mutation.counts), prob= eca.genome.length, replace=T)
        sum(eca.mutation.counts[res])*3.3*10^9/sum(eca.genome.length[res])
      })
      mutation.time.eca.lower <- c(mutation.time.eca.lower,quantile(bootstrapped.time, 0.025))
      mutation.time.eca.upper <- c(mutation.time.eca.upper,quantile(bootstrapped.time, 0.975))
      names(mutation.time.eca.lower)[length(mutation.time.eca.lower)] <- i
      names(mutation.time.eca.upper)[length(mutation.time.eca.upper)] <- i
    }else{
      mutation.time.eca.lower <- c(mutation.time.eca.lower,NA)
      mutation.time.eca.upper <- c(mutation.time.eca.upper,NA)
      names(mutation.time.eca.lower)[length(mutation.time.eca.lower)] <- i
      names(mutation.time.eca.upper)[length(mutation.time.eca.upper)] <- i
    }
  }
  
  
  ## now, test whether the individual fragments match with the eca peak: 
  if(length(eca.genome.length)>0){
    does.fragment.match <- apply(rbind(eca.mutation.counts, eca.genome.length), 2, function(x){
      if(length(eca.mutation.counts)==1){
        x[1] <- round(x[1])
      }
      if(x[1]==0 & sum(eca.mutation.counts)==0){
        return(1)
      }
      if(x[1] <= round(sum(eca.mutation.counts))* x[2]/sum(eca.genome.length)){
        test <- pnbinom(q = x[1], size = round(sum(eca.mutation.counts)), prob = round(sum(eca.mutation.counts))/(round(sum(eca.mutation.counts))*(1+x[2]/sum(eca.genome.length))))
      }else{
        test <- pnbinom(q = x[1], size = round(sum(eca.mutation.counts)), prob = round(sum(eca.mutation.counts))/(round(sum(eca.mutation.counts))*(1+x[2]/sum(eca.genome.length))), lower.tail = F)
      }
      test
    })
    does.fragment.match <- p.adjust(does.fragment.match)
    does.fragment.match <- sapply(does.fragment.match, function(x){
      if(x < 0.01){
        F
      }else{
        T
      }
    })
    ## Report fragments that conform and don't conform to a common origin
    remaining.aneuploidies.conforming.single.event[[i]] <- names(does.fragment.match[does.fragment.match])
    remaining.aneuploidies.not.conforming.single.event[[i]] <- names(does.fragment.match[!does.fragment.match])
    
    ## If some fragments do not match, check whether they request an even earlier time point and if so, determine
    if(length(does.fragment.match)>0 & any(!does.fragment.match)){
      fragments.before.eca <- which(eca.mutation.counts[!does.fragment.match]*3.3*10^9/eca.genome.length[!does.fragment.match] < mutation.time.eca[length(mutation.time.eca)])
      
      if(length(fragments.before.eca)>0){
        
        earliest.mutation.count <- eca.mutation.counts[!does.fragment.match][fragments.before.eca]
        earliest.genome.length <- eca.genome.length[!does.fragment.match][fragments.before.eca]
        
        earliest.mutation.time <- c(earliest.mutation.time, sum(eca.mutation.counts[!does.fragment.match][fragments.before.eca])*3.3*10^9/
                                      sum(eca.genome.length[!does.fragment.match][fragments.before.eca]))
        names(earliest.mutation.time)[length(earliest.mutation.time)] <- i
        
        ## bootstrap upper and lower limits of the mutation time
        if(length(earliest.mutation.count)>0){
          bootstrapped.time <- sapply(1:1000, function(x){
            res <- sample(x=1:length(earliest.mutation.count), size=length(earliest.mutation.count), prob= earliest.genome.length, replace=T)
            sum(earliest.mutation.count[res])*3.3*10^9/sum(earliest.genome.length[res])
          })
          earliest.mutation.time.lower <- c(earliest.mutation.time.lower,quantile(bootstrapped.time, 0.025))
          earliest.mutation.time.upper <- c(earliest.mutation.time.upper,quantile(bootstrapped.time, 0.975))
          names(earliest.mutation.time.lower)[length(earliest.mutation.time.lower)] <- i
          names(earliest.mutation.time.upper)[length(earliest.mutation.time.upper)] <- i
        }else{
          earliest.mutation.time.lower <- c(earliest.mutation.time.lower,NA)
          earliest.mutation.time.upper <- c(earliest.mutation.time.upper,NA)
          names(earliest.mutation.time.lower)[length(earliest.mutation.time.lower)] <- i
          names(earliest.mutation.time.upper)[length(earliest.mutation.time.upper)] <- i
        }
        
        
        does.fragment.match <- apply(rbind(earliest.mutation.count, earliest.genome.length), 2, function(x){
          if(length(earliest.mutation.count)==1){
            x[1] <- round(x[1])
          }
          if(x[1]==0 & sum(earliest.mutation.count)==0){
            return(1)
          }
          if(x[1] <= round(sum(earliest.mutation.count))* x[2]/sum(eca.genome.length)){
            test <- pnbinom(q = x[1], size = round(sum(earliest.mutation.count)), prob = round(sum(earliest.mutation.count))/(round(sum(earliest.mutation.count))*(1+x[2]/sum(earliest.genome.length))))
          }else{
            test <- pnbinom(q = x[1], size = round(sum(earliest.mutation.count)), prob = round(sum(earliest.mutation.count))/(round(sum(earliest.mutation.count))*(1+x[2]/sum(earliest.genome.length))), lower.tail = F)
          }
          test
        })
        does.fragment.match <- p.adjust(does.fragment.match)
        does.fragment.match <- sapply(does.fragment.match, function(x){
          if(x < 0.01){
            F
          }else{
            T
          }
        })
        ## Report fragments that conform and don't conform to a common origin
        aneuploidies.at.earliest.time[[i]] <- names(does.fragment.match[does.fragment.match])
        remaining.aneuploidies.not.conforming.earliest.event[[i]] <- names(does.fragment.match[!does.fragment.match])
      }
      
      
      
    }
    
    ## Finally, does an aneuploidy matching the MRCA alsow match with eca?
    if(length(eca.mutation.counts.at.mrca)>0){
      does.fragment.match <- apply(rbind(eca.mutation.counts.at.mrca, eca.genome.length.at.mrca), 2, function(x){
        if(length(eca.mutation.counts.at.mrca)==1){
          x[1] <- round(x[1])
        }
        if(x[1] < 1){
          x[1] <- round(x[1])
        }
        if(x[1]==0 & round(sum(eca.mutation.counts) + x[1])==0){return(T)}
        if(x[1] <= round((sum(eca.mutation.counts) + x[1]))* x[2]/(sum(eca.genome.length)+x[2]) ){
          test <- pnbinom(q = x[1], size = round(sum(eca.mutation.counts) + x[1]), prob = round(sum(eca.mutation.counts) + x[1])/(round(sum(eca.mutation.counts) + x[1])*(1+x[2]/(sum(eca.genome.length) + x[2]))))
        }else{
          test <- pnbinom(q = x[1], size = round(sum(eca.mutation.counts) + x[1]), prob = round(sum(eca.mutation.counts) + x[1])/(round(sum(eca.mutation.counts) + x[1])*(1+x[2]/(sum(eca.genome.length) + x[2]))), lower.tail = F)
        }
        test
      })
      does.fragment.match <- p.adjust(does.fragment.match)
      does.fragment.match <- sapply(does.fragment.match, function(x){
        if(x < 0.01){
          F
        }else{
          T
        }
      })
      aneuploidies.at.mrca.conforming.eca[[i]] <- names(does.fragment.match[does.fragment.match])
    }
    
  }
  
  names(mutation.time.eca)[length(mutation.time.eca)] <- i
}

### depending on your definition, define the ECA as the earliest event or not
mutation.time.eca[names(earliest.mutation.time)] <- earliest.mutation.time
mutation.time.eca.lower[names(earliest.mutation.time)] <- earliest.mutation.time.lower
mutation.time.eca.upper[names(earliest.mutation.time)] <- earliest.mutation.time.upper







