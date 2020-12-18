The script Mutation_density_quantification.R is a pipeline to estimate mutation densities at MRCA and ECA from tumor WG sequencing data. In the first part, mutations are classified as subclonal, non-amplified clonal or amplified clonal. In the second part, these estimates are being used to learn the densities at MRCA and ECA. 

The key output are the vectors ‘mutation.time.mrca’, ‘mutation.time.mrca.lower’,  ‘mutation.time.mrca.upper’,  ‘mutation.time.eca’,  ‘mutation.time.eca.upper’,  ‘mutation.time.eca.lower’, which contain the estimated mutation load per haploid genome at MRCA or ECA per tumor (and upper and lower 95% CI according to bootstrapping).
In addition, the following lists are returned 
- ‘aneuploidies.at.mrca’: chromosomal gains mapping to MRCA
- ‘remaining.aneuploidies.conforming.single.event’: chromosomal gains not mapping to MRCA, but to ECA
- ‘aneuploidies.at.mrca.conforming.eca’: chromosomal gains mapping to both MRCA and ECA
- ‘remaining.aneuploidies.not.conforming.single.event’: chromosomal gains neither mapping to MRCA nor to ECA
- ‘aneuploidies.at.earliest.time’: chromosomal gains mapping to a time point earlier than ECA

In order to execute the script, you need to set the working directory to the folder containing this script. In addition, you need to specify the path to the data and the IDs of the samples to analyze. The data should be organized in the same directory with subdirectories named by the tumor IDs. Each tumor subdirectory must contain two folders:
- ACEseq, containing the ACEseq output for this tumor (*comb_pro_extra*)
- SNVs, containing the .vcf file for this tumor. 

Example data are provided in the folder “Example_data”. Please set the data.directory to this folder and specify the tumors as 

tumors <- c("Example_tumor1_primary", "Example_tumor1_relapse", "Example_tumor2_primary",
            "Example_tumor3_primary", "Example_tumor3_relapse")
primary.of.tumor.pairs <- c("Example_tumor1_primary", "Example_tumor3_primary")
relapse.of.tumor.pairs <- c("Example_tumor1_relapse", "Example_tumor3_relapse")