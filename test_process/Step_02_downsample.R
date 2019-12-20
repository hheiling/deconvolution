## Using output from Step_01_PrepareBAM.R code, run downsample_prog.R code in IsoDeconvMM R package inst/ folder

#-------------------------------------------------------------------#
# Generate Downsampled Files:                                       #
#-------------------------------------------------------------------#

downsample_BAM<-function(infile,outlabel,directory,props.out, desct){
  cmd1 = sprintf("samtools view %s | wc -l > tmp_count.txt",infile)
  system(cmd1)
  # total_ct = scan("tmp_count.txt")
  
  for(j in props.out){
    # prob_val = (j/100)*(desct)/total_ct
    prob_val = (j/100)
    cmd2 = sprintf("java -jar /home/hheiling_unc_edu/picard.jar DownsampleSam INPUT=%s OUTPUT=%s%s_%s.bam PROBABILITY=%s VALIDATION_STRINGENCY=LENIENT RANDOM_SEED=null",infile,directory,outlabel,j,prob_val)
    cmd2
    system(cmd2)
  }
  
  # system("rm tmp_count.txt")
}

#---------------- APPLY DOWNSAMPLE_BAM -----------------------------#
# Set the directory where input files are located:
# inputDir = "/home/hheiling_unc_edu/deconvolution/test_process/test_materials"
# setwd(inputDir)

# Set Output Directory:
outDirec = "/home/hheiling_unc_edu/deconvolution/test_process/test_materials/"


# Insert Input Files:
inFiles = c("mm9_set1_sorted_by_name_uniq_filtered.bam",
            "mm9_set2_sorted_by_name_uniq_filtered.bam")

# Insert Output Labels:
outlabels = c("dsmm9_set1", "dsmm9_set2")

# Proportion Lists:
prop.vecs = list()
# prop.vecs[[1]] = seq(from = 10,to=100,by=10)
# prop.vecs[[2]] = seq(from = 10,to=100,by=10)
prop.vecs[[1]] = 50
prop.vecs[[2]] = 50

for(i in 1:length(inFiles)){
  cat(sprintf("start inFile: %s \n", inFiles[i]))
  downsample_BAM(infile = inFiles[i], outlabel = outlabels[i],
                 directory=outDirec, props.out=prop.vecs[[i]])
  cat(sprintf("end inFile: %s \n", inFiles[i]))
}

#----------------------------------------------------------#
# MERGING DOWNSAMPLED DATA                                 #
#----------------------------------------------------------#

merge_mat = matrix(c("dsmm9_set1_50.bam","dsmm9_set2_50.bam","mf_set1_50_set2_50.bam"),byrow=TRUE,ncol=3)

merge_func<-function(merge_mat,directory){
  for(i in 1:nrow(merge_mat)){
    cmd1 = sprintf("java -jar /home/hheiling_unc_edu/picard.jar MergeSamFiles I=%s I=%s O=%s/%s SORT_ORDER=queryname",merge_mat[i,1],merge_mat[i,2],directory,merge_mat[i,3])
    system(cmd1)
  }
}

#----------------------------------------------------------#
# APPLY THE FUNCTION                                       #
#----------------------------------------------------------#
directory = "/home/hheiling_unc_edu/deconvolution/test_process/test_materials/"

merge_func(merge_mat = merge_mat,directory = directory)
