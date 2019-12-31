## Using output from Step_01_PrepareBAM.R code, run downsample_prog.R code in IsoDeconvMM R package inst/ folder


# Needed arguments (see Step_0_Directions.Rmd for details):
# specifications

inFiles = specifications$inFiles
outlabels = specifications$outlabels
workingDirectory = specifications$workingDirectory
picard.jar.Dir = specifications$picard.jar.Dir
prop.mat = specifications$prop.mat
seeds = specifications$seeds
mix_outputs = mix_outputs


#-------------------------------------------------------------------#
# Generate Downsampled Files:                                       #
#-------------------------------------------------------------------#

downsample_BAM<-function(infile,outlabel,directory,props.out, rand_seeds){
  cmd1 = sprintf("samtools view %s | wc -l > tmp_count.txt",infile)
  system(cmd1)
  ds_names = character(length(props.out))
  
  for(j in props.out){
    prob_val = (j/100)
    seed = rand_seed[j]
    cmd2 = sprintf("java -jar %s/picard.jar DownsampleSam INPUT=%s OUTPUT=%s/%s_%s.bam PROBABILITY=%s VALIDATION_STRINGENCY=LENIENT RANDOM_SEED=%s",picard.jar.Dir, infile,directory,outlabel,j,prob_val,seed)
    cmd2
    system(cmd2)
    ds_names[j] = sprintf("%s_%s.bam", outlabel,j)
  }
  
  system("rm tmp_count.txt")
  
  return(ds_names)
}

#---------------- APPLY DOWNSAMPLE_BAM -----------------------------#

# Insert Output Labels:
# outlabels = c("dsmm9_set1", "dsmm9_set2")

ds_names_mat = matrix(NA, nrow = nrow(prop.mat), ncol = ncol(prop.mat))

for(i in 1:length(inFiles)){
  cat(sprintf("start inFile: %s \n", inFiles[i]))
  ds_names_mat[,i] = downsample_BAM(infile = inFiles[i], outlabel = outlabels[i],
                                    directory=workingDirectory, props.out=prop.mat[,i])
  cat(sprintf("end inFile: %s \n", inFiles[i]))
}

#----------------------------------------------------------#
# MERGING DOWNSAMPLED DATA                                 #
#----------------------------------------------------------#

# merge_mat = matrix(c("dsmm9_set1_50.bam","dsmm9_set2_50.bam","mf_set1_50_set2_50.bam"),byrow=TRUE,ncol=3)
merge_mat = cbind(ds_names_mat, mix_outputs)

merge_func<-function(merge_mat,directory){
  for(i in 1:nrow(merge_mat)){
    cmd1 = sprintf("java -jar %s/picard.jar MergeSamFiles I=%s I=%s O=%s/%s SORT_ORDER=queryname",picard.jar.Dir, merge_mat[i,1],merge_mat[i,2],directory,merge_mat[i,3])
    system(cmd1)
  }
}

#----------------------------------------------------------#
# APPLY THE FUNCTION                                       #
#----------------------------------------------------------#
directory = workingDirectory

merge_func(merge_mat = merge_mat,directory = directory)
