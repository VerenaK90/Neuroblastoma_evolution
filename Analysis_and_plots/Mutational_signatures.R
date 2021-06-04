##############################################################################################################################################
## Infer mutational signatures dominating SSNVs in neuroblastoma

##############################################################################################################################################
purity <- c()
ploidy <- c()
input.to.signatures <- data.frame(Sample = c(), chr = c(), pos = c(), ref = c(), alt = c(), subclonal = c())

for(i in tumors.80x){
  
  print(i)
  
  ## Read in ploidy-purity estimates:
  aceseq <- list.files(paste0(data.directory.80x, i, "/", cnv.directory, "/"), pattern="comb_pro_extra")[1]
  if(is.na(aceseq)){next}
  purity. <- Extract.purity.ploidy.from.ACEseq(aceseq)$purity
  ploidy. <- Extract.purity.ploidy.from.ACEseq(aceseq)$ploidy
  
  purity[i] <- purity.
  ploidy[i] <- ploidy.
  
  ## read in the mutation file
  files <- list.files(paste0(data.directory.80x, i, "/", snv.directory, "/"), pattern="somatic_snvs_conf_8_to_10")[1]
  
  mutations <- read.vcf(paste0(data.directory.80x, i, "/", snv.directory, "/", files))
  
  mutations$vcf$VAF <- Extract.info.from.vcf(mutations, info="VAF", type="snvs", mutationcaller="DKFZ")
  
  copy.number.info <- read.delim(file=paste0(data.directory.80x, i, "/", cnv.directory, "/", aceseq), sep="\t", stringsAsFactors = F)
  ## obtain the coverage ratios for the mutations of interest
  cnv.info.per.mutation <- Extract.copy.number.info.per.SSNV(mutations, copy.number.info)
  
  mutations <- as.data.frame(mutations$vcf)
  
  mutations$cov.ratio <- cnv.info.per.mutation$coverage.ratio
  
  ## classify mutations as clonal or subclonal based on a binomial distribution

  mutations$Subclonal <- apply(mutations, 1, function(x){
    CN <- round((as.numeric(x["cov.ratio"])*ploidy[i]*purity.+2*(1-purity[i])*(as.numeric(x["cov.ratio"])-1))/(purity[i]))
    if(is.na(x["cov.ratio"])|is.na(x["VAF"])|as.numeric(x["cov.ratio"])==0|round(as.numeric(x["cov.ratio"])*ploidy[i])==0){return(NA)}
    if(as.numeric(x["VAF"]) < qbinom(p = 0.05, size = 80, prob =  1*purity[i]/(purity[i]*CN + 2*(1-purity[i])))/80){
      "SC"}else{
        "C"
      }
  })
  
  
  
  tmp <- mutations
  tmp <- tmp[tmp$CHROM %in% c(1:22),]
  tmp <- data.frame(Sample = rep(i, nrow(tmp)), chr = tmp$CHROM, pos = tmp$POS, 
                    ref = tmp$REF, alt = tmp$ALT, subclonal = tmp$Subclonal)
  input.to.signatures <- rbind(input.to.signatures, tmp)
}

input.to.signatures <- input.to.signatures[!is.na(input.to.signatures$subclonal),]

input.to.signatures$chr <- paste0("chr", input.to.signatures$chr)

##############################################################################################################################################
# Convert to deconstructSigs input: all mutations
sigs.input <- mut.to.sigs.input(mut.ref = input.to.signatures, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg19)

# Determine the signatures contributing to the individual samples

# Plot output

signature.contribution <- matrix(0, nrow=30, ncol=nrow(sigs.input), dimnames = list(paste0("AC", 1:30), rownames(sigs.input)))

for(i in rownames(sigs.input)){
  out = whichSignatures(tumor.ref = sigs.input, 
                        signatures.ref = signatures.cosmic, 
                        sample.id = i, 
                        contexts.needed = TRUE,
                        tri.counts.method = 'genome',
                        signature.cutoff = 0.01)
  
  
  # Plot output
  plotSignatures(out, sub = i)
  signature.contribution[,i] <- unlist(out$weights)
  
}

##############################################################################################################################################
# Convert to deconstructSigs input: subclonal mutations
## now look at the mutational signature of subclonal mutations only 

sigs.input <- mut.to.sigs.input(mut.ref = input.to.signatures[input.to.signatures$subclonal=="SC",], 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg19)


# Determine the signatures contributing to the individual samples

# Plot output
signature.contribution.sc <- matrix(0, nrow=30, ncol=nrow(sigs.input), dimnames = list(paste0("AC", 1:30), rownames(sigs.input)))

for(i in rownames(sigs.input)){
  out = whichSignatures(tumor.ref = sigs.input, 
                        signatures.ref = signatures.cosmic, 
                        sample.id = i, 
                        contexts.needed = TRUE,
                        tri.counts.method = 'genome',
                        signature.cutoff = 0.01)
  
  signature.contribution.sc[,i] <- unlist(out$weights)
  
  # Plot output
  plotSignatures(out, sub = i)
}


##############################################################################################################################################
# Convert to deconstructSigs input: clonal mutations
sigs.input <- mut.to.sigs.input(mut.ref = input.to.signatures[input.to.signatures$subclonal!="SC",], 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = BSgenome.Hsapiens.UCSC.hg19)

signature.contribution.c <- matrix(0, nrow=30, ncol=nrow(sigs.input), dimnames = list(paste0("AC", 1:30), rownames(sigs.input)))

for(i in rownames(sigs.input)){
  out = whichSignatures(tumor.ref = sigs.input, 
                        signatures.ref = signatures.cosmic, 
                        sample.id = i, 
                        contexts.needed = TRUE,
                        tri.counts.method = 'genome',
                        signature.cutoff = 0.01)
  
  signature.contribution.c[,i] <- unlist(out$weights)
  
  
  # Plot output
  plotSignatures(out, sub = i)
}


##############################################################################################################################################


## do an overall picture, bootstrap the contribution to each signature

tmp <- data.frame(ALL=rep(0, 30),
                  SC = rep(0, 30),
                  C = rep(0, 30),
                  ALL.min = rep(0, 30),
                  ALL.max = rep(0, 30),
                  SC.min = rep(0, 30),
                  SC.max = rep(0, 30),
                  C.min = rep(0, 30),
                  C.max = rep(0, 30))

## overall picture:

tmp$ALL <- rowSums(signature.contribution)/sum(signature.contribution)
## error from bootstrapping
tmp.bootstrapped <- matrix(0, nrow(signature.contribution), 1000)
for(i in 1:1000){
  random.samples <- sample(1:ncol(signature.contribution), replace = T, size = ncol(signature.contribution))
  tmp.bootstrapped[,i] <-  rowSums(signature.contribution[,random.samples])/sum(signature.contribution[,random.samples])
}

tmp$ALL.min <- apply(tmp.bootstrapped, 1, quantile, p=0.025)
tmp$ALL.max <- apply(tmp.bootstrapped, 1, quantile, p=0.975)


tmp$C <- rowSums(signature.contribution.c)/sum(signature.contribution.c)
## error from bootstrapping
tmp.bootstrapped <- matrix(0, nrow(signature.contribution.c), 1000)
for(i in 1:1000){
  random.samples <- sample(1:ncol(signature.contribution.c), replace = T, size = ncol(signature.contribution.c))
  tmp.bootstrapped[,i] <-  rowSums(signature.contribution.c[,random.samples])/sum(signature.contribution.c[,random.samples])
}

tmp$C.min <- apply(tmp.bootstrapped, 1, quantile, p=0.025)
tmp$C.max <- apply(tmp.bootstrapped, 1, quantile, p=0.975)

tmp$SC <- rowSums(signature.contribution.sc)/sum(signature.contribution.sc)
## error from bootstrapping
tmp.bootstrapped <- matrix(0, nrow(signature.contribution.sc), 1000)
for(i in 1:1000){
  random.samples <- sample(1:ncol(signature.contribution.sc), replace = T, size = ncol(signature.contribution.sc))
  tmp.bootstrapped[,i] <-  rowSums(signature.contribution.sc[,random.samples])/sum(signature.contribution.sc[,random.samples])
}

tmp$SC.min <- apply(tmp.bootstrapped, 1, quantile, p=0.025)
tmp$SC.max <- apply(tmp.bootstrapped, 1, quantile, p=0.975)

