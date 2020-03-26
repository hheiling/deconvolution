# Run isoDu from isoform package

library(stringr)
library(isoform)

# Arrays 1-247
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

batch <- array_val
batch

# Only run 4 clusters at a time for each two-way cell type comparison
batch_seq = ((batch-1)*4 + 1):(batch*4)
batch_seq

# Home and scratch directories
header_nas = "/nas/longleaf/home/hheiling/Blueprint/"
header_pine = "/pine/scr/h/h/hheiling/Blueprint/"

# isoDetector output prefix
prefix_isoDetect = str_c(header_nas,"isoDetector_out")
# isoDu results prefix
prefix_isoDu = str_c(header_nas,"isoDu_out")
if(!dir.exists(prefix_isoDu)){dir.create(prefix_isoDu, recursive = T)}
# Human materials output prefix
prefix_mat = str_c(header_pine,"Blueprint_Materials")
# Fragment length prefix
prefix_fragLen = str_c(header_pine,"Fragment_Lengths")

# Datasets corresponding to cell types
CT = c("EGAD00001002671","EGAD00001002674","EGAD00001002675")
# Number of cell types
J = length(CT)

# Create output subdirectories
if(!dir.exists(str_c(prefix_isoDu, "CT1_vs_CT2", sep="/"))){
  dir.create(str_c(prefix_isoDu, "CT1_vs_CT2", sep="/"), recursive = T)
}
if(!dir.exists(str_c(prefix_isoDu, "CT1_vs_CT3", sep="/"))){
  dir.create(str_c(prefix_isoDu, "CT1_vs_CT3", sep="/"), recursive = T)
}
if(!dir.exists(str_c(prefix_isoDu, "CT2_vs_CT3", sep="/"))){
  dir.create(str_c(prefix_isoDu, "CT2_vs_CT3", sep="/"), recursive = T)
}


Routs = list()
samp_names = list()
for(ct in 1:length(CT)){
  Routs[[ct]] = list.files(path = str_c(prefix_isoDetect, CT[ct], sep="/"), pattern = ".RData", full.names = T)
  print(basename(Routs[[ct]]))
  samp_names[[ct]] = str_remove(basename(Routs[[ct]]), "_geneModel_knownIsoforms.RData")
  print(samp_names[[ct]])
}

Routs_list = list()
for(ct in 1:length(CT)){
  for(i in 1:length(Routs[[ct]])){
    # Load geneMod objects from isoDetector output
    load(Routs[[ct]][i])
    Routs_list[[i+(ct-1)*length(Routs[[ct]])]] = geneMod
  }
}

# Assume ordering of CT1, CT2, and CT3 corresponds to ordring of CT vector given above
tags =  c(rep("CT1",times=10),rep("CT2",times=10),rep("CT3",times=10))

# Get fragment size files
set.seed(6711)
fragFileList = list()
for(ct in 1:length(CT)){
  all_files = list.files(path = str_c(prefix_fragLen,CT[ct],sep="/"), 
                         pattern = "fragLen", full.names = T)
  samps = samp_names[[ct]]
  fragFiles = character(length(samps))
  for(s in 1:length(samps)){
    indic = sum(str_detect(all_files, samps[s]))
    if(indic == 1){
      fragFiles[s] = all_files[which(str_detect(all_files, samps[s]))]
    }else{
      fragFiles[s] = sample(all_files, size = 1)
    }
  }
  fragFileList[[ct]] = fragFiles
  print(sprintf("sample names for CT %s", CT[ct]))
  print(samps)
  print(sprintf("fragLen files for CT %s", CT[ct]))
  print(basename(fragFileList[[ct]]))
}

fragSizeFiles = c(fragFileList[[1]],fragFileList[[2]],fragFileList[[3]])

frag_list = readFragmentSizes(fragSizeFiles, lmax = 600)

# Read depth and restriction of clusters (those in isoDetector output)
readDepth = numeric(length(Routs_list))
clust_names = list()
for(i in 1:length(Routs_list)){
  geneMod = Routs_list[[i]]
  clust_names[[i]] = names(geneMod)
  cts = numeric(length(geneMod))
  for(g in 1:length(geneMod)){
    cts[g] = sum(geneMod[[g]]$y)
  }
  readDepth[i] = sum(cts)
}

restrict_clust = intersect(clust_names[[1]], clust_names[[2]])
for(i in 3:length(clust_names)){
  restrict_clust = intersect(restrict_clust, clust_names[[i]])
}

# load nTE data.frame
load(sprintf("%s/gencode.v15.nTE.RData", prefix_mat))
nTE[1:3,]

# Find genes and associated clusters to test for differential isoform usage
## Restrict to just a few chromosomes (not including X and Y sex chromosomes)
nTE_A = nTE[c(which(str_detect(nTE$clustID, "chr1_")),
              which(str_detect(nTE$clustID, "chr2_")),
              which(str_detect(nTE$clustID, "chr3_")),
              which(str_detect(nTE$clustID, "chr4_"))),]
## Restrict to clusters with number isoforms between 3 and 11
nTE_filtered = nTE_A[which(nTE_A$nT >= 3 & nTE_A$nT <= 11),]
dim(nTE_filtered)

# Combine isoDetector restriction with nTE filtering restrictions
restrict_clust = intersect(restrict_clust, nTE_filtered$clustID)

length(restrict_clust) # 986

# xData Covariates: Reference cell coding
xData = rep(c(0,1), each = 10)
print(xData)

if(any(batch_seq > length(restrict_clust))){
  if(all(batch_seq > length(restrict_clust))){
    stop("batch sequence not applicable to available clusters")
  }else{
    batch_seq = batch_seq[which(batch_seq <= length(restric_clust))]
  }
}

# Select small group of clusters to test at a time
clusts2use = restrict_clust[batch_seq]

# CT1 vs CT2
isoDu(tags = tags[1:20], gL = Routs_list[1:20], xData = xData, 
      outputFileName = sprintf("%s/CT1_vs_CT2/Batch_%s_isoDu_Usg.txt", prefix_isoDu, as.character(batch)),
      pdDistL = frag_list[1:20], g2test = clusts2use,
      readLen = 100, readDepth = readDepth[1:20], nResample = 200,
      method = "permutation", lmax = 600, duOnly = T)

# CT1 vs CT3
isoDu(tags = tags[c(1:10,21:30)], gL = Routs_list[c(1:10,21:30)], xData = xData, 
      outputFileName = sprintf("%s/CT1_vs_CT3/Batch_%s_isoDu_Usg.txt", prefix_isoDu, as.character(batch)),
      pdDistL = frag_list[c(1:10,21:30)], g2test = clusts2use,
      readLen = 100, readDepth = readDepth[c(1:10,21:30)], nResample = 200,
      method = "permutation", lmax = 600, duOnly = T)

# CT2 vs CT3
isoDu(tags = tags[11:30], gL = Routs_list[11:30], xData = xData, 
      outputFileName = sprintf("%s/CT2_vs_CT3/Batch_%s_isoDu_Usg.txt", prefix_isoDu, as.character(batch)),
      pdDistL = frag_list[11:30], g2test = clusts2use,
      readLen = 100, readDepth = readDepth[11:30], nResample = 200,
      method = "permutation", lmax = 600, duOnly = T)

print(gc())

# rm(list = ls())

#################################################################################################
