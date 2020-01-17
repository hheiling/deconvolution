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

# Combo Output Labels:
comboLabels = mix_labels


fragLengths<-function(Input_Files,outputLabels,comboList,comboLabels,useCombo){
  outLabels = paste(outputLabels,"_lengths.txt",sep="")
  for(i in 1:length(Input_Files)){
    cmd1 = sprintf("samtools view -f 65 %s | awk '{print ($8>=$4) ? $8-$4+%s : $4-$8+%s}' > %s",inputFiles[i], readLen, readLen, outLabels[i])
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
      cat(sprintf("comboList %i done \n", k))
    }
    
  }
}

fragLengths(Input_Files = inputFiles,outputLabels=outputLabels,comboList = comboList,comboLabels=comboLabels,useCombo=1)

# The End