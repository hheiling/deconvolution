
library(stringr)

prefix = "/pine/scr/h/h/hheiling/Blueprint"

# Download meta data
meta = readRDS(sprintf("%s/Blueprint_Materials/blueprint_meta_info.rds", prefix))
dim(meta)

# Identify outilers wrt paired vs single end reads
table(meta$LIBRARY_LAYOUT, meta$DATASET_ID) 

# From the above tabulation, remove the samples that had single end reads from the `EGAD00001002671`
# dataset and remove the samples that had paired end reads from the other two datasets. 

datasets = unique(meta$DATASET_ID)
datasets
lib_layout = c("SINGLE","SINGLE","PAIRED")

ls1 = list()
for(i in 1:length(datasets)){
  print(sprintf("Dataset: %s", datasets[i]))
  dat_sub = meta[which(meta$DATASET_ID == datasets[i] & meta$LIBRARY_LAYOUT == lib_layout[i]),]
  print(dim(dat_sub))
  ls1[[datasets[i]]] = dat_sub
}


set.seed(3407)

ls2 = list()

for(i in 1:length(ls1)){
  meta_sub = ls1[[i]]
  mix_rows = sample(1:nrow(meta_sub), size = round(nrow(meta_sub)/2))
  mix_sub = meta_sub[mix_rows,]
  fit_sub = meta_sub[-mix_rows,]
  
  print(names(ls1)[1])
  print(dim(fit_sub))
  print(dim(mix_sub))
  ls2[[names(ls1)[i]]] = list(mix_samps = mix_sub, fit_samps = fit_sub)
}


dirs_orig = list.dirs(str_c(prefix,"/Pure_ExonSetCts/"), full.names = T)
dirs_orig = dirs_orig[-1]
basename(dirs_orig)

for(i in 1:length(dirs_orig)){
  
  files_all = list.files(dirs_orig[i], full.names = T)
  print(sprintf("length of all files: %i", length(files_all)))
  meta_ls = ls2[[which(str_detect(names(ls2), basename(dirs_orig)[i]))]]
  
  # Identify files for fit
  fit_samps = meta_ls$fit_samps$EGAF
  print(sprintf("length fit_samps: %i", length(fit_samps)))
  fit_files = files_all[which(str_detect(basename(files_all),fit_samps[1]))]
  for(j in 2:length(fit_samps)){
    fit_files = union(fit_files, files_all[which(str_detect(basename(files_all),fit_samps[j]))])
  }
  print(sprintf("length fit_files: %i",length(fit_files)))
  print(head(basename(fit_files)))

  # Transfer copy of fit files to new directory
  fit_out_dir = str_c(prefix, "/Fit_Samples/", basename(dirs_orig)[i])

  for(f in 1:length(fit_files)){
    system(sprintf("cp %s %s", fit_files[f], fit_out_dir))
  }
  
  # Identify files for mixture creation
  mix_samps = meta_ls$mix_samps$EGAF
  mix_files = files_all[which(str_detect(basename(files_all),mix_samps[1]))]
  print(sprintf("length mix_samps: %i", length(mix_samps)))
  for(j in 2:length(mix_samps)){
    mix_files = union(mix_files, files_all[which(str_detect(basename(files_all),mix_samps[j]))])
  }
  print(sprintf("length mix_files: %i",length(mix_files)))
  
  # Transfer copy of mixture creation files to new directory
  mix_out_dir = str_c(prefix, "/MixCreation_Samples/", basename(dirs_orig)[i])
  
  for(f in 1:length(mix_files)){
    system(sprintf("cp %s %s", mix_files[f], mix_out_dir))
  }
}

# Check results

for(dir1 in c("Fit_Samples","MixCreation_Samples")){
  for(dir2 in datasets){
    files = list.files(path = str_c(prefix, "/", dir1, "/", dir2))
    print(sprintf("Number of samples in %s/%s directory: %i", dir1, dir2, length(files)))
  }
}


##############################################################################################