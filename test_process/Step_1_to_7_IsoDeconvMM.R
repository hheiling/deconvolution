## Running of actual IsoDeconvMM R package

# Required Libraries:
library(gtools)
# library(IsoDeconvMM)
# source("https://bioconductor.org/biocLite.R")
# biocLite("cummeRbund")
library(cummeRbund)
# Had to downgrade to RSQ-lite v 1.1-2 --> ?
library(alabama)
library(stringr)

prequel = "/home/hheiling_unc_edu/IsoDeconvMM/R/"
to_source = c("isoDeconv_geneModel_revised.R",
              "geneModel_multcell_edit.R",
              "loadData_edit.R",
              "pdist_gen.R",
              "Production_Functions_MixedSamp.R",
              "Production_Functions_PureSamp.R",
              "rem_clust.R",
              "EffectiveLength_Total.R")
source_files = str_c(prequel, to_source)

for(i in 1:length(to_source)){
  source(file = source_files[i])
}



# Inputs
file = "set1_50_set2_50" # Should be one of the comboLabels provided in FragLengths_Function_mixture.R

sys_statement1 = sprintf("countData=c(\"mf_%s_counts.txt\",
                         \"mm9_set1_counts.txt\",
                         \"mm9_set2_counts.txt\")",file)
eval(parse(text=sys_statement1))
labels = c("mix","set1_ref1", "set2_ref1")
cellTypes = c("mix","set1","set2")
fragSizeFile = sprintf("%s_fraglens.txt",file)
bedFile = "Mus_musculus.NCBIM37.67.nonoverlap.exon.bed"
knownIsoforms = "Mus_musculus.NCBIM37.67.nonoverlap.exon.knownIsoforms.RData"
readLen=76
eLenMin=1

lmax=600

# folder = "/home/hheiling_unc_edu/deconvolution/test_process/test_materials"
total_cts = numeric(length(countData))
for(i in 1:length(countData)){
  countsi = read.table(countData[i], as.is = T)
  counts_col = countsi[,1]
  total_cts[i] = sum(counts_col)
}

# Step 1

files = file

final_geneMod = list()

for(j in 1:length(files)){
  # Call dev_compiled_geneMod function
  fin_geneMod = dev_compiled_geneMod(countData=countData,labels = labels,total_cts = total_cts, 
                                     cellTypes=cellTypes, bedFile=bedFile,knownIsoforms=knownIsoforms,
                                     fragSizeFile=fragSizeFile,readLen=readLen,lmax=lmax,
                                     eLenMin=eLenMin)
  final_geneMod[[j]] = fin_geneMod
}

# Warning message:
#   In fin_geneMod["Sample_Info"] <- list(info = info_mat, tclust_tot = length(fin_geneMod),  :
#                                           number of items to replace is not a multiple of replacement length

save(final_geneMod, file = "Step1_final_geneMod.RData")

# Step 2

#-----------------------------------------------------------------------------#
# ESTABLISHING CLUSTERS WITH HIGHEST LEVELS OF DISCRIMINATORY CAPABILITY      #
#-----------------------------------------------------------------------------#

#------------------ Identify Highly Discriminatory Clusters -------------------#

# Pick genes from chromosome 18 with 3 or more isoforms (transcripts)

# siggenes object from "Exploring RData Objects.R"
# Columns: geneId, nE, nT, clustID
# Folder: deconvolution/test_process/test_materials/
load("chr18_siggenes.RData")

finalOut2 = siggenes

# Step 3

analy_genes = finalOut2$geneId

significant_geneMod = list()

for(j in 1:length(final_geneMod)){
  
  fin_geneMod = final_geneMod[[j]]
  
  indices2chk = which(names(fin_geneMod)!="Sample_Info")
  indices_tmp = NULL
  indices=NULL
  # indices_tmp = rep(0,length(geneMod)) # What is this geneMod object?
  # Idea: should be fin_geneMod instead of geneMod
  indices_tmp = rep(0,length(fin_geneMod))
  for(i in indices2chk){
    infodf = fin_geneMod[[i]]$info
    genesi = unique(infodf$gene)
    genesi = unique(unlist(strsplit(x=genesi,split = ":")))
    if(any(genesi %in% analy_genes)){indices_tmp[i]=1}
  }
  indices = which(indices_tmp==1)
  
  sig_geneMod = fin_geneMod[indices]
  sig_geneMod["Sample_Info"] = fin_geneMod["Sample_Info"]
  
  sig_geneMod = rem_clust(geneMod = sig_geneMod,co = 5,min_ind = 0)
  
  significant_geneMod[[j]] = sig_geneMod
}

save(significant_geneMod, file = "Step3_significant_geneMod.RData")

# Step 4

#-------------------------------------------------------------------#
# EDIT TO GROUP CELL TYPES                                          #
#-------------------------------------------------------------------#

modified_sig_geneMod = list()

for(f in 1:length(significant_geneMod)){
  
  sig_geneMod = significant_geneMod[[f]]
  
  info_mat = sig_geneMod[["Sample_Info"]]
  cellTypes = unique(info_mat$Cell_Type)
  
  ctList = list()
  
  for(j in 1:length(cellTypes)){
    idx = which(info_mat$Cell_Type==cellTypes[j])
    ctList[[cellTypes[j]]] = list(samps = info_mat$Label[idx], tots = info_mat$Total[idx])
  }
  
  idx2consider = which(names(sig_geneMod)!="Sample_Info")
  for(k in idx2consider){
    for(l in 1:length(cellTypes)){
      samps2use = ctList[[l]]$samps
      tots      = ctList[[l]]$tots
      
      y_vecs  = paste("sig_geneMod[[k]]$y",samps2use,sep = "_")
      y_vecsc = paste(y_vecs,collapse = ",")
      nExon = eval(parse(text=sprintf("length(%s)",y_vecs[1])))
      textcmd = sprintf("matrix(c(%s),nrow=nExon,ncol=length(samps2use))",y_vecsc)
      expMat  = eval(parse(text=textcmd))
      
      totmg   = tots-colSums(expMat)
      expMat2 = rbind(totmg,expMat)
      
      if(cellTypes[l]!="mix"){
        sig_geneMod[[k]][[cellTypes[l]]] = list(cellType=cellTypes[l],rds_exons=expMat2)
      } else {
        sig_geneMod[[k]][[cellTypes[l]]] = list(cellType=cellTypes[l],rds_exons_t=expMat2)
      }
      
    }
  }
  
  modified_sig_geneMod[[f]] = sig_geneMod
  
}

save(modified_sig_geneMod, file = "Step4_modified_sig_geneMod.RData")

# Step 5

#-----------------------------------------------------------#
# CALL Pure Sample                                          #
#-----------------------------------------------------------#

## Need further investiation here / in above steps ##

cellTypes = c("set1","set2")

pure_est = list()

for(j in 1:length(modified_sig_geneMod)){
  
  sig_geneMod = modified_sig_geneMod[[j]]
  
  sim.out = sig_geneMod[which(names(sig_geneMod)!="Sample_Info")]
  
  # Clusters with single isoforms:
  # EXCLUDE THEM FOR THE MOMENT!.
  dim_mat = matrix(0,nrow=length(sim.out),ncol=2)
  excl_clust = c()
  excl_clust2 = c()
  
  # for(i in 1:(length(sim.out))){ # See explanation of change below
  for(i in 1:(length(sim.out)-1)){
    dim_mat[i,] = dim(sim.out[[i]][["X"]])
    # dim_mat[i,] = dim(sim.out[[i]]["X"][[1]])
    
    if(all(dim_mat[i,]==c(1,1))){
      excl_clust = c(excl_clust,i)
    }
    if(dim_mat[i,2] == 1){
      excl_clust2 = c(excl_clust2,i)
    }
  }
  
  # Note: Exclude last sim.out entry due to explanation below
  excl_clust2 = c(excl_clust2, length(sim.out))
  
  excl_clust_union = union(excl_clust,excl_clust2)
  if(length(excl_clust_union)>0){
    sim.new = sim.out[-c(excl_clust_union)]
  } else {
    sim.new = sim.out
  }
  
  # Optimize the Pure Sample Functions:
  tmp.data = Pure.apply.fun(data.list = sim.new, cellTypes = cellTypes, corr_co = 1)
  
  pure_est[[j]] = tmp.data
  
}


# Error in `[.data.frame`(sim.out[[i]], "X") : undefined columns selected
  # Error understood due to comment below

# Note: Very last entry in sim.out (sim.out[[182]] in this example) was not the usual "list"
# object that the other sim.out elements had.
# For now, removing this object, but should investigate later

# > tmp.data = Pure.apply.fun(data.list = sim.new, cellTypes = cellTypes, corr_co = 1)
# Error in base::colSums(x, na.rm = na.rm, dims = dims, ...) : 
#   'x' must be an array of at least two dimensions

# Solution: Found that sim.new[[i]][cellTypes]$setk$rds_exons had one column for all elements,
# so colSums function would not work. Corrected code to accommodate matrix with only 1 column.

save(pure_est, file = "Step5_pure_est.RData")

# Step 6

IsoDeconv_Output = list()

for(i in 1:length(pure_est)){
  
  tmp.data = pure_est[[i]]
  
  #--------------------------------------------------------#
  # Establish input break ups                              #
  #--------------------------------------------------------#
  # Cell-Types:
  cellTypes = c("set1","set2")
  
  # Data Set Necessities:
  clust.start = 1
  clust.end = length(tmp.data)
  by.value = 15
  
  start.pts = seq(from = 1,to = clust.end,by = by.value)
  end.pts = c((start.pts[-1]-1),clust.end)
  
  cluster_output = list()
  for(m in 1:length(start.pts)){
    start.pt = start.pts[m]
    end.pt = end.pts[m]
    # Call Revised_Sim_MixCode_SI.R code
    curr.clust.opt = tmp.data[c(start.pt:end.pt)]
    curr.clust.out = STG.Update_Cluster.All(all_data=curr.clust.opt, cellTypes = cellTypes,
                                            optimType="nlminb", simple.Init=FALSE, initPts=c(0.5))
    
    cluster_output[[m]] = curr.clust.out
  }
  
  IsoDeconv_Output[[i]] = cluster_output
}

# Simple Init Not Performed! Full Fit for each start point!
#   Error in cdata[["alpha.est"]][, cellTypes] : subscript out of bounds

# Fixed problem by correcting cellTypes to be c("set1","set2")

save(IsoDeconv_Output, file = "Step6_IsoDeconv_Output.RData")

# Step 7

#-----------------------------------------------------------------------#
# Compile Files                                                         #
#-----------------------------------------------------------------------#

Final_Compiled_Output = list()

for(j in 1:length(IsoDeconv_Output)){
  
  comp.out= NULL
  curr.clust.out = NULL
  
  #---- Set up new pattern ----#
  est.chunks = IsoDeconv_Output[[j]]
  
  message("File ", j)
  
  #---- Populate Output Dataset ----#
  comp.out = list()
  r = 1
  
  for(i in 1:length(est.chunks)){
    curr.clust.out = est.chunks[[i]]
    nl = length(curr.clust.out)
    for(m in 1:nl){
      comp.out[[r]]=curr.clust.out[[m]]
      r = r+1
    }
  }
  
  Final_Compiled_Output[[j]] = comp.out 
  
}

save(Final_Compiled_Output, file = "Step7_FinalCompiled_Output.RData")
