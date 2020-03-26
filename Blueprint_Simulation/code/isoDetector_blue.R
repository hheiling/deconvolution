# Running isoDetector for each pure reference sample

library(stringr)
library(isoform)

# Arrays 1-30
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
batch <- array_val

header_nas = "/nas/longleaf/home/hheiling/Blueprint/"
header_pine = "/pine/scr/h/h/hheiling/Blueprint/"

# Chose cell type (dataset)
if(batch <= 10){
  CT = "EGAD00001002671"
  k = 1
  j = batch
}else if(batch >= 21){
  CT = "EGAD00001002675"
  k = 3
  j = batch - 20 # want j from 1 to 10 
}else{
  CT = "EGAD00001002674"
  k = 2
  j = batch - 10
}

# Set seed for random selection of pure samples
set.seed(8280)
seeds = sample(1000:9999, size = 3, replace = F)

# Simulated output prefix
prefix_out = str_c(header_nas,"isoDetector_out/",CT)
if(!dir.exists(prefix_out)){dir.create(prefix_out, recursive = T)}
# Fragment length prefix
prefix_fragLen = str_c(header_pine,"Fragment_Lengths/",CT)
# Materials prefix
prefix_mat = str_c(header_pine,"Blueprint_Materials")
# Pure cell type files prefix (to be used for algorithm fit)
prefix_pure = str_c(header_pine,"Fit_Samples/",CT)


# Find pure reference samples for CT (cell type)
pure_all = list.files(prefix_pure, pattern = "trec_exon_set.txt", full.names = T)
length(pure_all)

# Find fragment length files
fragSizeFiles = list.files(path = prefix_fragLen, pattern = "fragLen.txt", full.names = T)
length(fragSizeFiles)

# Select sampling of pure files and fragment length files
set.seed(seeds[k])
pure_files = sample(pure_all, size = 10, replace = F)
print(basename(pure_files))

pure_select = pure_files[j] # Select specific file to run isoDetector on
print(basename(pure_select))

samp_name = unlist(str_split(basename(pure_select), "_"))[1]
print(samp_name)

if(sum(str_detect(basename(fragSizeFiles), samp_name)) == 1){
  fragSizeFile = fragSizeFiles[which(str_detect(basename(fragSizeFiles), samp_name))]
  print(basename(fragSizeFile))
}else{
  fragSizeFile = sample(fragSizeFiles, size = 1)
  print(basename(fragSizeFile))
}

# Find BED and knownIsoforms objects
bedFile = sprintf("%s/gencode.v15.nonoverlap.exon.bed", prefix_mat)
knownIsoforms = sprintf("%s/gencode.v15.nonoverlap.exon.knownIsoforms.RData", prefix_mat)

# Specify name of output file
output_file = str_c(samp_name, "_geneModel_knownIsoforms.RData")

print(output_file)

isoDetector(pure_select, bedFile, fragSizeFile, readLen = 100, # readLen from step1_get_counts.Rout
            sprintf("%s/%s", prefix_out, output_file), 
            knownIsoforms=knownIsoforms)

