
workDir   = "/fh/scratch/delete90/sun_w/plittle/CS_eQTL/s5_EGA"
dataset   = "EGAD00001002675"
sam_names = list.files(file.path(workDir, dataset), pattern="EGAF")
length(sam_names)
sam_names[1:5]

# ------------------------------------------------------------------------
# read in sample information
# ------------------------------------------------------------------------

meta = readRDS("data/blueprint_meta_info.rds")
dim(meta)
meta[1:2,]

table(sam_names %in% meta$EGAF)

# ------------------------------------------------------------------------
# generiate bash file
# ------------------------------------------------------------------------

sh_file = sprintf("step1_get_counts_%s.sh", dataset)
cat("", file=sh_file)

for(sam_name in sam_names){
  cmd = "sbatch -c 6 --mem=22528 --partition=largenode --wrap="
  cmd = paste0(cmd, "\"R CMD BATCH --no-save --no-restore")
  cmd = sprintf("%s \\\"--args dataset='%s' sam_name='%s'\\\"", cmd, dataset, sam_name)
  cmd = sprintf("%s step1_get_counts.R step1_get_counts_%s_%s.Rout\"\n", 
                cmd, dataset, sam_name)
  cat(cmd, file=sh_file, append=TRUE)
}

gc()

sessionInfo()
q(save="no")
