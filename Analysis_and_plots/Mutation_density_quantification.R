##############################################################################################################################################
## Quantify the densitiy of amplified and non-amplified clonal mutations
##############################################################################################################################################
## Estimate the number of false positive clonal SNVs from the two sample pairs

source(paste0(custom.script.directory, "Estimate_sampling_bias_from_tumor_pairs.R"))
load(paste0(rdata.directory, "Assess_false_positives_due_to_sampling.RData"))

#source(paste0(custom.script.directory, "Adjust_purity.R"))
load(paste0(rdata.directory, "/Purity_ploidy.RData"))
##############################################################################################################################################
## For each of the tumors, we want to know the number of clonal mutations that lie on one or multiple copies of a gained chromosome.
## Clonal mutations on several copies were likely acquired prior to tumor initiation, given that the gain is clonal.
## Since the number is only informative if corrected for the respective genome fraction, we need to report this as well.


## Store the clonal mutations present on any number of alleles and store this information for each chromosome separately

clonal.mutations.all.tumors <- list()
clonal.cpg.mutations.all.tumors <- list()

## Iterate through the tumors, iterate through the autosomes and iterate through the different B allele frequencies

## First loop: tumors
for(i in c(tumors.discovery, tumors.validation)){

  print(i)

  if(i %in% tumors.discovery){
    data.directory <- data.directory.discovery
  }else{
    data.directory <- data.directory.validation
  }

  ## Find ACEseq file:
  aceseq <- list.files(paste0(data.directory, "/", i, "/", cnv.directory, "/"), pattern="comb_pro_extra")[1]
  if(is.na(aceseq)){next}
  aceseq <- paste0(data.directory, "/", i, "/", cnv.directory, "/", aceseq)

  ## Find mutation file
  mutations <- list.files(paste0(data.directory, "/", i, "/", snv.directory, "/"), pattern="somatic_snvs_conf_8_to_10")[1]
  mutations <- paste0(data.directory, "/", i, "/", snv.directory, "/", mutations)

  clonal.mutations.all.tumors[[i]] <- count.clonal.mutations(aceseq, mutations, chromosomes = c(1:22), purity.=purity[i], ploidy.=ploidy[i])
  clonal.cpg.mutations.all.tumors[[i]] <- count.clonal.mutations(aceseq, mutations, chromosomes = c(1:22), CpG_TpG = T, purity.=purity[i], ploidy.=ploidy[i])

}

B.allele.indicator <- clonal.mutations.all.tumors[[1]]$B.allele.indicator
copy.number.indicator <- clonal.mutations.all.tumors[[1]]$ copy.number.indicator

save(clonal.mutations.all.tumors, clonal.cpg.mutations.all.tumors,
     file=paste0(rdata.directory, "Clonal_mutations_different_CNs.RData"))

load(paste0(rdata.directory, "Clonal_mutations_different_CNs.RData"))


##############################################################################################################################################
## Compute the mutational density at each genomic fragment

## We model mutation accumulation as a Poisson process. The number of mutations acquired on a piece of DNA depends on the per-base mutation rate (at this time) and the length of the piece
## using a chisquare-approximation, we obtain confidence intervals for the mutation time.

## store results of all tumors in a list
mutation.time <- list()

for(i in c(tumors.discovery, tumors.validation)){

  mutation.time[[i]] <- data.frame()

  ## iterate through the chromosomes
  for(j in 1:22){

    ## Get normalized mutation counts per copy and segment on this chromosome
    tmp.mut.count <- Mutation.time.converter(clonal.mutations.all.tumors[[i]]$clonal.mutation.matrix[,j])

    tmp.genome.length <- clonal.mutations.all.tumors[[i]]$segment.length.matrix[,j]
    ## Restrict analysis to fragments > 10^7 bp
    tmp.mut.count <- tmp.mut.count[which(tmp.genome.length>10^7)]
    tmp.genome.length <- tmp.genome.length[tmp.genome.length>10^7]

    if(length(tmp.mut.count)==0){next}

    ## Convert to mutations per haploid genome and store output

    mutation.time[[i]] = rbind(mutation.time[[i]],
                          data.frame(Mean =  tmp.mut.count*3.3*10^9/tmp.genome.length,
                                     Min = 0.5*qchisq(0.025, tmp.mut.count*2)/tmp.genome.length*3.3*10^9,
                                     Max = 0.5*qchisq(0.975, (tmp.mut.count*2+2))/tmp.genome.length*3.3*10^9,
                                     Segment = paste("chr", j, names(tmp.mut.count), sep="_")))

  }

  mutation.time[[i]]$Segment <- factor(mutation.time[[i]]$Segment, levels = mutation.time[[i]]$Segment[order(mutation.time[[i]]$Mean)])

  ## Visualize the mutation density per segment
  p <- ggplot(mutation.time[[i]],
              aes(x=Segment, y=Mean, ymin=Min, ymax=Max)) + geom_col() + geom_errorbar(aes(x=Segment, ymin=Min, ymax=Max)) +
    scale_y_continuous(limits=c(0, max(mutation.time[[i]]$Max)), name = "# Mutations per haploid genome") +
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
earliest.mutation.time <- c()
mutation.time.eca <- c()

mrca.eca <- list()

for(i in c(tumors.discovery, tumors.validation)){
  print(i)

  ## build input matrix

  mrca.eca[[i]] <- MRCA.ECA.quantification(clonal.mutations.all.tumors[[i]]$clonal.mutation.matrix,
                                      clonal.mutations.all.tumors[[i]]$segment.length.matrix)

  mutation.time.mrca <- rbind(mutation.time.mrca, data.frame(Mean = mrca.eca[[i]]$mutation.time.mrca,
                                                             Min = mrca.eca[[i]]$mutation.time.mrca.lower,
                                                             Max = mrca.eca[[i]]$mutation.time.mrca.upper,
                                                             Sample = i))

  if(length(mrca.eca[[i]]$earliest.mutation.time)>0){
    earliest.mutation.time <- rbind(earliest.mutation.time, data.frame(Mean = mrca.eca[[i]]$earliest.mutation.time,
                                                               Min = mrca.eca[[i]]$earliest.mutation.time.lower,
                                                               Max = mrca.eca[[i]]$earliest.mutation.time.upper,
                                                               Sample = i))
  }

  mutation.time.eca <- rbind(mutation.time.eca, data.frame(Mean = mrca.eca[[i]]$mutation.time.eca,
                                                             Min = mrca.eca[[i]]$mutation.time.eca.lower,
                                                             Max = mrca.eca[[i]]$mutation.time.eca.upper,
                                                             Sample = i))

}

rownames(mutation.time.mrca) <- mutation.time.mrca$Sample
rownames(mutation.time.eca) <- mutation.time.eca$Sample
rownames(earliest.mutation.time) <- earliest.mutation.time$Sample

save(mutation.time.mrca, mutation.time.eca, earliest.mutation.time, mrca.eca, mutation.time,
     file=paste0(rdata.directory, "MRCA_timing.RData"))


median(mutation.time.mrca$Mean)/3.3/10^3
median(mutation.time.eca$Mean, na.rm=T)/3.3/10^3
quantile(mutation.time.mrca$Mean)/3.3/10^3
quantile(mutation.time.mrca$Mean[!is.na(mutation.time.eca$Mean)])/3.3/10^3
quantile(mutation.time.eca$Mean, na.rm=T)/3.3/10^3

mutation.time.eca[earliest.mutation.time$Sample,]$Mean <- earliest.mutation.time$Mean
mutation.time.eca[earliest.mutation.time$Sample,]$Min <- earliest.mutation.time$Min
mutation.time.eca[earliest.mutation.time$Sample,]$Max <- earliest.mutation.time$Max
## report this value
quantile(mutation.time.eca$Mean, na.rm=T)/3.3/10^3

## 0.008, 0.022, 0.048

##############################################################################################################################################
## Alternative check: do segment CIs overlap with MRCA?
evidence.for.eca.from.individual.ci <- c()

for(i in c(tumors.discovery, tumors.validation)){
  fragments.timing.mrca <- mutation.time[[i]]$Min[sapply(mutation.time[[i]]$Segment, function(x){
    x <- strsplit(as.character(x), split="_")[[1]][4]
    if(is.na(x) || x=="I"){T}else{
      F
    }
  })]
  evidence.for.eca.from.individual.ci[i] <- any(mutation.time[[i]]$Max < min(fragments.timing.mrca))
}
names(evidence.for.eca.from.individual.ci) <- c(tumors.discovery, tumors.validation)

intersect(names(evidence.for.eca.from.individual.ci==T), rownames(mutation.time.eca)[!is.na(mutation.time.eca$Mean)])
setdiff(names(evidence.for.eca.from.individual.ci==T), rownames(mutation.time.eca)[!is.na(mutation.time.eca$Mean)])
setdiff(rownames(mutation.time.eca)[!is.na(mutation.time.eca$Mean)], names(evidence.for.eca.from.individual.ci==T))
## would increase number of tumors with ECA

##############################################################################################################################################
### Extract for each tumor the VAF distribution for each ploidy state. Exclude sex chromosomes

vafs.all.tumors <- list()
genome.size.all.tumors <- list()

load(paste0(rdata.directory, "Purity_ploidy.RData"))

purities.all.tumors <- purity
ploidies.all.tumors <- ploidy
for(i in c(tumors.discovery, tumors.validation)){

  print(i)

  if(i %in% tumors.discovery){
    data.directory <- data.directory.discovery
  }else{
    data.directory <- data.directory.validation
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

  vafs.this.tumor <- list()
  genome.size.this.tumor <- list()
  p <- list()

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


       p[[length(p)+1]] <-  ggplot(data.frame(VAF=readcounts.[,2]/rowSums(readcounts.)), aes(x=VAF)) + geom_histogram(binwidth = 0.01) +
          geom_vline(xintercept = prob.clonal, col="firebrick", linetype=2) + geom_vline(xintercept=monosomic.prob.clonal, col="firebrick", linetype=2)+
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(name="# Mutations") +
          scale_x_continuous(limits=c(0,1)) + ggtitle(paste0("Copy number = ", k, "B allele = ", l))



      }
    }

    pdf(paste0(data.directory, i, "/VAF_histograms.pdf"), width=5, height=5)
    print(ggarrange(plotlist = p, nrow=2, ncol=2))
    dev.off()


  vafs.all.tumors[[i]] <- vafs.this.tumor
  genome.size.all.tumors[[i]] <- genome.size.this.tumor
  save(vafs.all.tumors, genome.size.all.tumors, file=paste0(rdata.directory, "Vafs_all_tumors.RData"))

}

save(vafs.all.tumors, genome.size.all.tumors, file=paste0(rdata.directory, "Vafs_all_tumors.RData"))



