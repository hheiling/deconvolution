## Step 0.3 described in IsoDeconvMM Pipline.Rmd file in inst/ folder of IsoDeconvMM R package

#--------------------------------------------------------------------------------------
# Set Parameters
#--------------------------------------------------------------------------------------

# Needed arguments (see Step_0_Directions.Rmd for details):
# workingDirectory, bedFile, sim (logical), and possibly mixture_files if sim = T

# Input directory: has the "sorted_by_name_uniq_filtered.bam" files
inputDir = workingDirectory # Output folder from Step_01 - Step_03

#------------------------------------------------------------------------------------------------
# Loop across all replicates and files:
#        For loop assumes that all files have same bed File.
#------------------------------------------------------------------------------------------------


## Counts for "_sorted_by_name_uniq_filtered.bam" files created in Step 01

cmd  = sprintf("ls %s/*_sorted_by_name_uniq_filtered.bam", workingDirectory)
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

## Simulated Mixture Samples

if(sim == TRUE){
  
  cmd  = sprintf("ls %s/%s", workingDirectory, mixture_files)
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
  
}

