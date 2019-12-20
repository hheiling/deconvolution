## Step 0.5 in IsoDeconvMM Pipeline.Rmd file in inst/ folder of IsoDeconvMM R package
## See also: FragLengths_Function_Mixture.R 

#-------------------------------------------------------------#
# Fragment Length Files: Generation                           #
#-------------------------------------------------------------#

# Specify Directory where read files are kept:
# setwd("deconvolution/test_process/test_materials")

# BAM Files where sequenced reads are located:
inputFiles = c("mf_set1_50_set2_50.bam",
               "mm9_set1_sorted_by_name_uniq_filtered.bam",
               "mm9_set2_sorted_by_name_uniq_filtered.bam")

# Associated Labels for output files:
outputLabels = c("set1_50_set2_50",
                 "set1_r1",
                 "set2_r1")

# If any files are to be combined, list them in separate units:
comboList = list()
comboList[[1]] = c("set1_r1_lengths.txt",
                   "set2_r1_lengths.txt",
                   "set1_50_set2_50_lengths.txt")

# Combo Output Labels:
comboLabels = c("set1_50_set2_50")


fragLengths<-function(Input_Files,outputLabels,comboList,comboLabels,useCombo){
  outLabels = paste(outputLabels,"_lengths.txt",sep="")
  for(i in 1:length(Input_Files)){
    cmd1 = sprintf("samtools view -f 65 %s | awk '{print ($8>=$4) ? $8-$4+76 : $4-$8+76}' > %s",inputFiles[i],outLabels[i])
    system(cmd1)
    cat(sprintf("Input_File %s samtools command done \n", Input_Files[i]))
  }
  
  if(missing(comboList) && useCombo==0){
    for(j in 1:length(outLabels)){
      cmd2_a = sprintf("cat %s | sort -n | uniq -c > %s_fraglens.txt",outLabels[j],outputLabels[j])
      system(cmd2_a)
    }
  } else if(useCombo==1 && missing(comboList)){
    stop("If you wish to combine files, you must list which files are to be combined!")
    
  } else { # comboList specified, and useCombo = 1
    
    ftc = unlist(lapply(X = comboList,FUN = function(x) {return(paste(x,collapse=" "))}))
    
    for(k in 1:length(comboList)){
      cat(sprintf("start comboList %i \n", k))
      cmd2_b = sprintf("cat %s | sort -n | uniq -c > %s_fraglens.txt",ftc[k],comboLabels[k])
      system(cmd2_b)
      cat(sprintf("comboList %i done \n"), k)
    }
    
  }
}

fragLengths(Input_Files = inputFiles,outputLabels=outputLabels,comboList = comboList,comboLabels=comboLabels,useCombo=1)

# The End