##############################################################################################################################################
## Script to extract functional mutations and likely driver mutations from the detection tumor cohort

##############################################################################################################################################
## Load settings

source("Settings.R")

##############################################################################################################################################
## For each of the tumors, extract coding mutations; structural variations, high-level amplifications and deletions

functional.mutations <- data.frame(CHROM=c(), POS=c(), ID=c(), REF=c(), ALT=c(), QUAL=c(), FIlTER=c(), INFO=c(), FORMAT=c(), FORMAT_INFO=c(),
                                   ANNOVAR_FUNCTION=c(), GENE=c(), EXONIC_CLASSIFICATION=c(), READS_REF=c(), READS_ALT=c(), 
                                   CoverageRatio = c(), BAF = c(), Purity = c(), Ploidy = c(), SAMPLE=c(), AA_change=c())

## store to generate a table later on
translocations.for.output <- c()
amplifications.for.output <- c()
deletions.for.output <- c()
wb <- createWorkbook()

translocations <- data.frame(gene1=c(), gene2=c(), Sample=c(), svtype=c())
amplifications <- data.frame(gene=c(), Sample=c())
deletions <- data.frame(gene=c(), Sample=c())

## iterate through the tumors
for(i in tumors.discovery){
  
  print(i)
  
  ## Read in combined estimates of ploidy and purity:
  aceseq <- list.files(paste0(data.directory.discovery,  i, "/", cnv.directory, "/"), pattern="comb_pro_extra")[1]
  if(is.na(aceseq)){next}
  purity. <- Extract.purity.ploidy.from.ACEseq(aceseq)$purity
  ploidy. <- Extract.purity.ploidy.from.ACEseq(aceseq)$ploidy
  
   
  ## Read in the mutation file
  files <- list.files(paste0(data.directory.discovery,  i, "/", snv.directory, "/"), pattern="somatic_snvs_conf_8_to_10")[1]
  if(is.na(files)){
    print(i)
    next
  }
  
  mutations <- read.vcf(paste0(data.directory.discovery,  i, "/", snv.directory, "/", files))
  mutations$vcf <- mutations$vcf[mutations$vcf$ANNOVAR_FUNCTION %in% c("exonic", "upstream", "splicing"),]
  colnames(mutations$vcf)[10] <- "FORMAT_INFO"
  mutations$vcf <- mutations$vcf[,c(1:10, 16, 17, 18, 19)]
  
  mutations$vcf$AA_change <- Extract.info.from.vcf(mutations, info="AA_change")

  
  counts <- Extract.info.from.vcf(mutations, info="readcounts")
  
  
  copy.number.info. <- read.delim(file=paste0(data.directory.discovery,  i, "/", cnv.directory, "/", aceseq), sep="\t", stringsAsFactors = F)
  ## obtain the coverage ratios/BAF/genotype/TCN for the mutations of interest
  cnv.info.per.mutation <- Extract.copy.number.info.per.SSNV(mutations, copy.number.info.)
  ## obtain the coverage ratios at mutated loci
  coverage.ratios. <- cnv.info.per.mutation$coverage.ratio
  bafs. <- cnv.info.per.mutation$baf
  A <- cnv.info.per.mutation$A
  copy.number. <- cnv.info.per.mutation$tcn
  
  mutations <- mutations$vcf
  mutations$READS_REF <- counts[,1]
  mutations$READS_ALT <- counts[,2]
  
  mutations$CoverageRatio <- coverage.ratios.
  mutations$BAF <- bafs.
  mutations$Purity <- purity.
  mutations$Ploidy <- ploidy.
  mutations$TCN <- copy.number.
  mutations$A <- A
  mutations$SAMPLE <- i
  mutations$ANNOVAR_TRANSCRIPTS <- NULL
  
  
  functional.mutations <- rbind(functional.mutations, mutations)
  
  
  ## same for indels
  
  ## read in the mutation file
  files <- list.files(paste0(data.directory.discovery,  i, "/", indel.directory, "/"), pattern="somatic_indels_conf_8_to_10")[1]
  if(is.na(files)){
    print(i)
    next
  }
  
  mutations <- read.vcf(paste0(data.directory.discovery,  i, "/", indel.directory, "/", files))
  mutations$vcf <- mutations$vcf[,c(1:9, 11, 17, 18, 19, 20)]
  colnames(mutations$vcf)[10] <- "FORMAT_INFO"
  mutations$vcf$AA_change <- Extract.info.from.vcf(mutations, info="AA_change", type = "indel")
  
  
  if(nrow(mutations$vcf)==0){next}
  
  counts <- Extract.info.from.vcf(mutations, info="readcounts", type = "indel")
  
  copy.number.info. <- read.delim(file=paste0(data.directory.discovery,  i, "/", cnv.directory, "/", aceseq), sep="\t", stringsAsFactors = F)
  ## obtain the coverage ratios for the mutations of interest
  cnv.info.per.mutation <- Extract.copy.number.info.per.SSNV(mutations, copy.number.info.)
  ## obtain the coverage ratios at mutated loci
  coverage.ratios. <- cnv.info.per.mutation$coverage.ratio
  bafs. <- cnv.info.per.mutation$baf
  A <- cnv.info.per.mutation$A
  copy.number. <- cnv.info.per.mutation$tcn
  
  mutations <- mutations$vcf
  mutations$READS_REF <- counts[,1]
  mutations$READS_ALT <- counts[,2]
  mutations$CoverageRatio <- coverage.ratios.
  mutations$BAF <- bafs.
  mutations$Purity <- purity.
  mutations$Ploidy <- ploidy.
  mutations$TCN <- copy.number.
  mutations$A <- A
  mutations$SAMPLE <- i
  mutations$ANNOVAR_TRANSCRIPTS <- NULL
  
  
  functional.mutations <- rbind(functional.mutations, mutations)
  
  ## look up structural variants 
  
  ## read in the mutation file
  files <- list.files(paste0(data.directory.discovery,  i, "/", sv.directory, "/"), pattern="filtered_somatic_minEventScore3.tsv")[1]
  if(is.na(files)){
    print(i)
    next
  }
  
  mutations <- read.vcf(paste0(data.directory.discovery,  i, "/", sv.directory, "/", files))
  mutations <- mutations$vcf
  mutations$genepos1 <- sapply(mutations$gene1, function(x){strsplit(x, split="_")[[1]][2]})
  mutations$genepos2 <- sapply(mutations$gene2, function(x){strsplit(x, split="_")[[1]][2]})
  
  mutations$gene1 <- sapply(mutations$gene1, function(x){strsplit(x, split="_")[[1]][1]})
  mutations$gene2 <- sapply(mutations$gene2, function(x){strsplit(x, split="_")[[1]][1]})
  
  ## Take translocations between 2 genes or within 1 gene, but not just intronic deletions/amplifications
  mutations <- mutations[mutations$gene1 != mutations$gene2 | mutations$genepos1 != mutations$genepos2,]
  ## Require an event score of >= 5
  mutations <- mutations[(mutations$gene1 %in% driver.genes | mutations$gene2 %in% driver.genes) &
                           as.numeric(mutations$eventScore) >= 5,,drop=F]
  if(nrow(mutations)>0){
    mutations$Sample <- i
    translocations <- rbind(translocations, mutations[,c("gene1", "gene2", "Sample", "svtype")])
    translocations.for.output <- rbind(translocations.for.output, mutations)
  }
  
  ### Look for high-level amplifications and deletions
  
  amp <- unlist(apply(driver.genes.with.genetic.pos, 1, function(x){
    tmp <- copy.number.info.[paste0("chr",copy.number.info.$X.chromosome)==x[2] & ((as.numeric(copy.number.info.$start) > as.numeric(x[3]) & as.numeric(x[4]) > as.numeric(copy.number.info.$start))|
                                                                                     (as.numeric(copy.number.info.$end) < as.numeric(x[4]) & as.numeric(x[3]) < as.numeric(copy.number.info.$end))|
                                                                                     (as.numeric(copy.number.info.$end) > as.numeric(x[4]) & as.numeric(x[3]) > as.numeric(copy.number.info.$start)) ), ,drop=F]
    
    if(nrow(tmp)>0){
      ## Require at least 10 copies
      if(mean(as.numeric(tmp$tcnMean))>=10 ){
        return(x[1])
      }else{
        return(c())
      }
    }else{
      return(c())
    }
  }))
  

  if(length(amp)>0){
    amplifications <- rbind(amplifications, data.frame(gene=amp, Sample=i))
    
    amp.for.output <- do.call("rbind", apply(driver.genes.with.genetic.pos, 1, function(x){
      tmp <- copy.number.info.[paste0("chr",copy.number.info.$X.chromosome)==x[2] & ((as.numeric(copy.number.info.$start) > as.numeric(x[3]) & as.numeric(x[4]) > as.numeric(copy.number.info.$start))|
                                                                                       (as.numeric(copy.number.info.$end) < as.numeric(x[4]) & as.numeric(x[3]) < as.numeric(copy.number.info.$end))|
                                                                                       (as.numeric(copy.number.info.$end) > as.numeric(x[4]) & as.numeric(x[3]) > as.numeric(copy.number.info.$start)) ), ,drop=F]
      
      if(nrow(tmp)>0){
        ## Require at least 10 copies
        if(mean(as.numeric(tmp$tcnMean))>=10 ){
          tmp$Gene <- x[1]
          tmp <- tmp[tmp$tcnMean>=10,]
          return(tmp)
        }else{
          return()
        }
      }else{
        return()
      }
    }))
    
    amp.for.output$Sample <- i
    amplifications.for.output <- rbind(amplifications.for.output, amp.for.output)
  }
  
  ## Homozygous deletions; require <0.9 copy numbers
  del <- unlist(apply(driver.genes.with.genetic.pos, 1, function(x){
    tmp <- copy.number.info.[paste0("chr",copy.number.info.$X.chromosome)==x[2] & ((as.numeric(copy.number.info.$start) > as.numeric(x[3]) & as.numeric(x[4]) > as.numeric(copy.number.info.$start))|
                                                                                     (as.numeric(copy.number.info.$end) < as.numeric(x[4]) & as.numeric(x[3]) < as.numeric(copy.number.info.$end))|
                                                                                     (as.numeric(copy.number.info.$end) > as.numeric(x[4]) & as.numeric(x[3]) > as.numeric(copy.number.info.$start)) ), ,drop=F]
    
    if(nrow(tmp)>0){
      if(mean(as.numeric(tmp$tcnMean))<=0.9){
        return(x[1])
      }else{
        return(c())
      }
    }else{
      return(c())
    }
  }))
  
 
  
  if(length(del)>0){
    deletions <- rbind(deletions, data.frame(gene=del, Sample=i))
    
    del.for.output <- do.call("rbind", apply(driver.genes.with.genetic.pos, 1, function(x){
      tmp <- copy.number.info.[paste0("chr",copy.number.info.$X.chromosome)==x[2] & ((as.numeric(copy.number.info.$start) > as.numeric(x[3]) & as.numeric(x[4]) > as.numeric(copy.number.info.$start))|
                                                                                       (as.numeric(copy.number.info.$end) < as.numeric(x[4]) & as.numeric(x[3]) < as.numeric(copy.number.info.$end))|
                                                                                       (as.numeric(copy.number.info.$end) > as.numeric(x[4]) & as.numeric(x[3]) > as.numeric(copy.number.info.$start)) ), ,drop=F]
      
      if(nrow(tmp)>0){
        if(mean(as.numeric(tmp$tcnMean))<=0.9){
          tmp$Gene <- x[1]
          tmp <- tmp[tmp$tcnMean<=0.9,]
          
          return(tmp)
        }else{
          return()
        }
      }else{
        return()
      }
    }))
    del.for.output$Sample <- i
    deletions.for.output <- rbind(deletions.for.output, del.for.output)
  }
  
  
}

translocations <- unique(translocations)
translocations$SV <- "SV"

## filter non-promoter TERT mutations
functional.mutations <- functional.mutations[-which(functional.mutations$GENE=="TERT" &
                                                      functional.mutations$SAMPLE=="XI003_19142"),]

## keep splice-site mutations, nonsynonymous SNVs and indels

to.keep <- functional.mutations$ANNOVAR_FUNCTION=="splicing" | 
(  functional.mutations$ANNOVAR_FUNCTION=="exonic" &
  functional.mutations$EXONIC_CLASSIFICATION != "synonymous SNV" )|
  functional.mutations$EXONIC_CLASSIFICATION =="ncRNA_exonic" |
  (functional.mutations$ANNOVAR_FUNCTION=="upstream" &
     functional.mutations$GENE=="TERT" & functional.mutations$POS %in% c(1295228, 1295250))


filtered.functional.mutations <- functional.mutations[to.keep,]

## annotate whether the mutation is a known cancer mutation

filtered.functional.mutations$Known_driver <- apply(filtered.functional.mutations, 1, function(x){
  tmp <- known.driver.positions[known.driver.positions$gene==x["GENE"],,drop=F]
  if(nrow(tmp)==0){
    return(F)
  }else if(any(tmp$amino_acid==paste("p", x["AA_change"], sep="."))){
    return(T)
  }else{
    return(F)
  }
})



## Ignore structural variants in ATRX as we have the more detailed info here already from Hartlieb et al., 2020
translocations <- unique(translocations)
translocations$SV <- "SV"
translocations <- translocations[translocations$gene1!="ATRX" & gene2!="ATRX",]
translocations.for.output <- translocations.for.output[translocations.for.output$gene1!="ATRX" & gene2!="ATRX",]

## Ignore MYCN as we have the more detailed info here already from Hartlieb et al., 2020 
setdiff(amplifications[amplifications$gene=="MYCN", "Sample"], names(telomere.classification.discovery[telomere.classification.discovery=="MNA"]))
setdiff( names(telomere.classification.discovery[telomere.classification.discovery=="MNA"]),amplifications[amplifications$gene=="MYCN", "Sample"])

amplifications <- amplifications[amplifications$gene!="MYCN",]
amplifications.for.output <- amplifications.for.output[amplifications.for.output$Gene!="MYCN",]

## take out ATRX, BCORL1, MED12; they lie on the x-chromosome
deletions <- deletions[!deletions$gene %in% c("ATRX", "MED12", "BCORL1"),]
deletions.for.output <- deletions.for.output[!deletions.for.output$Gene %in% c("ATRX", "MED12", "BCORL1"),]


## simplify by adding splice sites to Exonic classification
filtered.functional.mutations$EXONIC_CLASSIFICATION[is.na(filtered.functional.mutations$EXONIC_CLASSIFICATION)] <- filtered.functional.mutations$ANNOVAR_FUNCTION[is.na(filtered.functional.mutations$EXONIC_CLASSIFICATION)]
filtered.functional.mutations$Known_driver <- replace(filtered.functional.mutations$Known_driver, filtered.functional.mutations$Known_driver==T, "Known_driver")
filtered.functional.mutations$Known_driver <- replace(filtered.functional.mutations$Known_driver, filtered.functional.mutations$Known_driver==F, "")

## manually add as Known driver: frameshift deletions/stopgain in intogen candidates and ALK mutation F1174S and mutations in ATRX
filtered.functional.mutations[filtered.functional.mutations$GENE %in% driver.genes &
                                filtered.functional.mutations$EXONIC_CLASSIFICATION %in% c("stopgain", "stoploss", "frameshift deletion", "frameshift insertion"),]$Known_driver <- "Known_driver"

filtered.functional.mutations[filtered.functional.mutations$AA_change=="F1174S",]$Known_driver <- "Known_driver"
filtered.functional.mutations[filtered.functional.mutations$GENE=="ATRX",]$Known_driver <- "Known_driver"

##############################################################################################################################################
## transform the mutation info into a mutation matrix, where rows correspond to genes/chromosomes and columns to samples
mat <- matrix("", nrow=length(unique(driver.genes)), ncol=length(tumors.discovery),
              dimnames=list(c(unique(driver.genes)), tumors.discovery))

translocations <- translocations[,c("gene1", "gene2", "Sample", "SV")]
transolcations <- unique(translocations)

for(i in 1:nrow(mat)){
    mat[i,as.character(filtered.functional.mutations[filtered.functional.mutations$GENE==rownames(mat)[i],]$SAMPLE)] <- filtered.functional.mutations[filtered.functional.mutations$GENE==rownames(mat)[i],]$EXONIC_CLASSIFICATION
    
    goi <- rownames(mat)[i]
    mat[i,as.character(deletions[deletions$gene==goi,]$Sample)] <- paste(mat[i,as.character(deletions[deletions$gene==goi,]$Sample)], "DEL", sep=";")
    mat[i,as.character(amplifications[amplifications$gene==goi,]$Sample)] <- paste(mat[goi,as.character(amplifications[amplifications$gene==goi,]$Sample)], "AMP", sep=";")
    mat[i,as.character(translocations[translocations$gene1==goi | translocations$gene2==goi,]$Sample)] <- paste(mat[goi,as.character(translocations[translocations$gene1==i | translocations$gene2==goi,]$Sample)], translocations[translocations$gene1==goi | translocations$gene2==goi,]$SV, sep=";")
   # mat[i,as.character(filtered.functional.mutations[filtered.functional.mutations$GENE==rownames(mat)[i],]$SAMPLE)] <- paste(mat[i,as.character(filtered.functional.mutations[filtered.functional.mutations$GENE==rownames(mat)[i],]$SAMPLE)], filtered.functional.mutations[filtered.functional.mutations$GENE==rownames(mat)[i],]$Known_driver, sep=";")
}


## Manual adjustments:
## We don't consider deletions in MYCN as drivers
mat["MYCN",][mat["MYCN",]==";DEL"] <- ""
mat <- mat[rowSums(mat!="")>=0,]
## Take information from Hartlieb et al., 2020 for ATRX
mat["ATRX",] <- sample.information.discovery[colnames(mat), "ATRX"]
mat[is.na(mat)] <- ""
mat[mat=="wt"] <- ""
## Don't consider whole chromosome loss of chromosome X
mat[mat=="whole chromosome loss"] <- ""
mat[mat=="whole chromosome loss; stopgain"] <- "stopgain"
## Take information from Hartlieb et al., 2020 for MNA/TERT
mat["MYCN",sample.information.discovery[colnames(mat),"Telomere.maintenance.mechanism"]=="MNA"] <- paste(mat["MYCN",sample.information.discovery[colnames(mat),"Telomere.maintenance.mechanism"]=="MNA"], "AMP", sep=";")
mat["TERT",sample.information.discovery[colnames(mat),"Telomere.maintenance.mechanism"]=="TERT"] <- paste(mat["TERT",sample.information.discovery[colnames(mat),"Telomere.maintenance.mechanism"]=="TERT"], "SV", sep=";")
## adjust heterogeneous cases (they were annnotated as "HET" in the excel file but have SVs in TERT and amplifications in MYCN)
mat["TERT","NBE57"] <- "SV"
mat["MYCN","NBE57"] <- "AMP"
mat["MYCN","NBE19"] <- "AMP"
mat["TERT","NBE85"] <- "SV"
