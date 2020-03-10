## Step 0.5 in IsoDeconvMM Pipeline.Rmd file in inst/ folder of IsoDeconvMM R package
## See also: FragLengths_Function_Mixture.R 

#-------------------------------------------------------------#
# Fragment Length Files: Generation                           #
#-------------------------------------------------------------#

# Needed arguments (see Step_0_Directions.Rmd for details):
# mixture_files, pure_files, mix_labels, pure_labels, readLen

# BAM Files where sequenced reads are located:
inputFiles = c(mixture_files, pure_files)

# Associated Labels for output files:
outputLabels = c(mix_labels, pure_labels)

# If any files are to be combined, list them in separate units:
comboList = list()
for(i in 1:length(mixture_files)){
  comboList[[i]] = c(str_c(pure_labels, "_lengths.txt"), str_c(mix_labels[i], "_lengths.txt"))
}

Input_Files  = "/Users/wsun/research/data/blueprint/bams_processed/EGAF00001329931_sorted_by_name_uniq_filtered.bam"
outputLabels = "/Users/wsun/research/data/blueprint/bams_processed/EGAF00001329931"
readLen = 100

# Combo Output Labels:
comboLabels = mix_labels

outLabels = paste(outputLabels,"_lengths.txt",sep="")

for(i in 1:length(Input_Files)){
  cmd1 = sprintf("samtools view -f 65 %s | awk '{print ($8>=$4) ? $8-$4+%s : $4-$8+%s}' > %s",
                 Input_Files[i], readLen, readLen, outLabels[i])
  system(cmd1)
  cat(sprintf("Input_File %s samtools command done \n", Input_Files[i]))
}

for(j in 1:length(outLabels)){
  cmd2_a = sprintf("cat %s | sort -n | uniq -c > %s_fraglens.txt",outLabels[j], outputLabels[j])
  system(cmd2_a)
}


# The End