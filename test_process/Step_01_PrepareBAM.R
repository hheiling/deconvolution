## Step 0.1 Described in "IsoDeconvMM Pipeline.Rmd" in IsoDeconvMM R package inst/ folder

library(asSeq)

# inputDir = "deconvolution/test_process/"
inputFolder = "/home/hheiling_unc_edu/deconvolution/test_process"
outputFolder = "/home/hheiling_unc_edu/deconvolution/test_process/test_materials"

# Check cat file | less to see if appropriate output given (to exit, hit q)

# Keep inputDir and outputFolder equal to the same directory:
# outputFolder = "/home/hheiling_unc_edu/deconvolution/test_process"
# When I did this, the original "mm9_set1_sorted_by_name.bam" files were replaced or renamed with the 
# new "_sorted_by_name_unique_filtered.bam" prefix, but at least this gave good output

#------------------------------------------------------------------
# Generating BAM file List            
#------------------------------------------------------------------

#set input directory
# setwd(inputDir)

#Generate initial list of files: (NOTE: assume .bam.bai files also included)
init_list = list.files(pattern=".bam")

#Generate list of .bai files
int_list = list.files(pattern=".bai")

#Final List (excluding .bai files):
BAM_list = setdiff(init_list,int_list)

# -----------------------------------------------------------------
# check the bam files
# -----------------------------------------------------------------

#Checks length of BAM files to ensure all has run properly:
length(BAM_list)

#Displays BAM list as another check for errors:
BAM_list

bam2use = BAM_list

# nchar("_sorted_by_name") = 15

#Loop across all BAM files in the folder:
for(i in 1:length(BAM_list)){
  
  bami = bam2use[i]
  
  #Generates a name that can be used in the files by stripping 
  #off the ".bam" extension: (4+15 = 19 characters)
  sami = substr(bam2use[i],start=1,stop=nchar(bam2use[i])-4) 
  
  # ----------------------------------------------------------
  # counting
  # ----------------------------------------------------------
  ctF  = sprintf("%s/count_%s.txt", inputFolder,sami)
  cmd1 = sprintf("samtools view %s | wc -l >> %s \n", bam2use[i], ctF)
  system(cmd1)
  
  # ----------------------------------------------------------
  # sorting
  # ----------------------------------------------------------
  cmd2 = sprintf("samtools sort -n %s -o %s/%s", bam2use[i], inputFolder, sami)
  system(cmd2)
  bamF = sprintf("%s/%s.bam", outputFolder, sami)
  
  # ----------------------------------------------------------
  # getUnique and filtering
  # ----------------------------------------------------------
  prepareBAM(bamF, sprintf("%s/%s", inputFolder, sami), sortIt=FALSE)
  
  system(sprintf("rm %s", bamF))
  
  # ----------------------------------------------------------
  # counting again
  # ----------------------------------------------------------
  cmd3   = sprintf("samtools view %s/%s_sorted_by_name_uniq_filtered.bam | wc -l >> %s\n", inputFolder, sami, ctF)
  system(cmd3)
  
}

# ----------------------------------------------------------
# sorting by position
# ----------------------------------------------------------
#Generate initial list of files: (NOTE: assume .bam.bai files also included)
to_sort_list = list.files(pattern="sorted_by_name_uniq_filtered.bam")

#Checks length of BAM files to ensure all has run properly:
length(to_sort_list)

#Displays BAM list as another check for errors:
to_sort_list

for(i in 1:length(to_sort_list)){
  cmd4 = sprintf("samtools sort %s -o %s/%s", to_sort_list[i], outputFolder, to_sort_list[i])
  system(cmd4)
}

# Experimentation:

outputFolder = "/home/hheiling_unc_edu/Step_01_test"
i=1
  
  bami = bam2use[i]

  #Generates a name that can be used in the files by stripping
  #off the ".bam" extension: (4+15 = 19 characters)
  sami = substr(bam2use[i],start=1,stop=nchar(bam2use[i])-4)
  # 
  # # ----------------------------------------------------------
  # # counting
  # # ----------------------------------------------------------
  # ctF  = sprintf("%s/count_%s.txt", outputFolder,sami)
  # cmd1 = sprintf("samtools view %s | wc -l >> %s\n", bam2use[i], ctF)
  # system(cmd1)
  
  system("samtools view mm9_set1_sorted_by_name.bam |wc -l >> /home/hheiling_unc_edu/Step_01_test/counts_mm9_set1_sorted_by_name.txt")
  system("samtools view mm9_set2_sorted_by_name.bam |wc -l >> /home/hheiling_unc_edu/Step_01_test/counts_mm9_set2_sorted_by_name.txt")

  # ----------------------------------------------------------
  # sorting
  # ----------------------------------------------------------
  # cmd2 = sprintf("samtools sort -n %s -o %s/%s", bam2use[i], outputFolder, sami)
  # system(cmd2)
  # bamF = sprintf("%s/%s.bam", outputFolder, sami)
  
  system("samtools sort -n mm9_set1_sorted_by_name.bam -o /home/hheiling_unc_edu/Step_01_test/mm9_set1_sorted_by_name")
  system("samtools sort -n mm9_set2_sorted_by_name.bam -o /home/hheiling_unc_edu/Step_01_test/mm9_set2_sorted_by_name")
  
  # ----------------------------------------------------------
  # getUnique and filtering
  # ----------------------------------------------------------
  prepareBAM(bamF, sprintf("%s/%s", outputFolder, sami), sortIt=FALSE)
  
  prepareBAM("/home/hheiling_unc_edu/Step_01_test/mm9_set1_sorted_by_name.bam", 
             "/home/hheiling_unc_edu/Step_01_test/mm9_set1_sorted_by_name", 
             sortIt = FALSE)
  
  prepareBAM("/home/hheiling_unc_edu/Step_01_test/mm9_set2_sorted_by_name.bam", 
             "/home/hheiling_unc_edu/Step_01_test/mm9_set2_sorted_by_name", 
             sortIt = FALSE)
  
  # system(sprintf("rm %s", bamF))
  
  # ----------------------------------------------------------
  # counting again
  # ----------------------------------------------------------
  # cmd3   = sprintf("samtools view %s/%s_sorted_by_name_uniq_filtered.bam | wc -l >> %s\n", outputFolder, sami, ctF)
  # system(cmd3)
  
  system("samtools view mm9_set1_sorted_by_name_uniq_filtered.bam |wc -l >> /home/hheiling_unc_edu/Step_01_test/counts_mm9_set1_sorted_by_name.txt")
  system("samtools view mm9_set2_sorted_by_name_uniq_filtered.bam |wc -l >> /home/hheiling_unc_edu/Step_01_test/counts_mm9_set2_sorted_by_name.txt")
  

  # ----------------------------------------------------------
  # sorting again - by position
  # ----------------------------------------------------------
  system("samtools sort mm9_set1_sorted_by_name_uniq_filtered.bam -o /home/hheiling_unc_edu/Step_01_test/Step_01_output/mm9_set1_sorted_by_name_uniq_filtered.bam")
  system("samtools sort mm9_set2_sorted_by_name_uniq_filtered.bam -o /home/hheiling_unc_edu/Step_01_test/Step_01_output/mm9_set2_sorted_by_name_uniq_filtered.bam")
  