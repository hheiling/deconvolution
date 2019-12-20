## Step 0.3 described in IsoDeconvMM Pipline.Rmd file in inst/ folder of IsoDeconvMM R package

#--------------------------------------------------------------------------------------
# Set Parameters
#--------------------------------------------------------------------------------------

# Input directory: has the "sorted_by_name_uniq_filtered.bam" files
inputDir = "/home/hheiling_unc_edu/deconvolution/test_process/test_materials" # Output folder from Step_01 - Step_03

# BED file (once identified and downloaded)
bedFile = "Mus_musculus.NCBIM37.67.nonoverlap.exon.bed"

#--------------------------------------------------------------------------------------
# ISOFORM Software Library Access
#--------------------------------------------------------------------------------------
library(isoform) 

#------------------------------------------------------------------------------------------------
# Set Working Directory
#------------------------------------------------------------------------------------------------
# setwd(inputDir)

#------------------------------------------------------------------------------------------------
# Loop across all replicates and files:
#        For loop assumes that all files within a class (ECC1, HMEC, GM12878) have same bed File.
#------------------------------------------------------------------------------------------------


## "Pure" cell type samples

cmd  = "ls *_sorted_by_name_uniq_filtered.bam"
ffs  = system(cmd, intern=TRUE)
length(ffs)
head(ffs)
sams = gsub("_sorted_by_name_uniq_filtered.bam", "", ffs)

for(i in 1:length(ffs)){
  
  sam1 = sams[i]
  cat(i, sam1, date(), "\n")
  
  bamFile = ffs[i]
  outFile = sprintf("%s_counts.txt", sam1)
  
  countReads(bamFile, bedFile, outFile)
}

## Mixture samples

cmd  = "ls mf_set1_50_set2_50.bam"
ffs  = system(cmd, intern=TRUE)
length(ffs)
head(ffs)
sams = gsub(".bam", "", ffs)

for(i in 1:length(ffs)){
  
  sam1 = sams[i]
  cat(i, sam1, date(), "\n")
  
  bamFile = ffs[i]
  outFile = sprintf("%s_counts.txt", sam1)
  
  countReads(bamFile, bedFile, outFile)
}