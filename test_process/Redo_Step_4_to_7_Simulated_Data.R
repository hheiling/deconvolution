# Re-testing code using simulated data

# Load step4_sim_data() list object
# load("test_process/test_materials/Step4_Simulated_Additions.RData")

library(alabama)

prequel = "C:/Users/hheiling/Documents/GitHub/IsoDeconvMM/R/"
to_source = c("Production_Functions_MixedSamp.R",
              "Production_Functions_PureSamp.R")
source_files = str_c(prequel, to_source)

for(i in 1:length(to_source)){
  source(file = source_files[i])
}

modified_sig_geneMod2 = list()
modified_sig_geneMod2[[1]] = step4_sim_data

# Step 5

#-----------------------------------------------------------#
# CALL Pure Sample                                          #
#-----------------------------------------------------------#

## Need further investiation here / in above steps ##

cellTypes = c("set1","set2")

pure_est = list()

for(j in 1:length(modified_sig_geneMod2)){
  
  sig_geneMod = modified_sig_geneMod2[[j]]
  
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

save(Final_Compiled_Output, file = "test_process/test_materials/Step7_Redo_Final_Compiled_Output.RData")

#-------------------------------------------------------------------#
# Exploration of Output                                             #
#-------------------------------------------------------------------#

# Load Final_Compiled_Output object from above output
load("test_process/test_materials/Step7_Redo_Final_Compiled_Output.RData")

explore = Final_Compiled_Output[[1]]
length(explore)
explore_1 = explore[[1]]
length(explore_1)
names(explore_1$mix)
explore_1$mix$p.est

# Summary of p.est values

p.est.1 = numeric(length(explore))
p.est.2 = numeric(length(explore))

for(i in 1:length(explore)){
  p.est.1[i] = explore[[i]]$mix$p.est[,1]
  p.est.2[i] = explore[[i]]$mix$p.est[,2]
}

final_p.est = data.frame(cellType1 = p.est.1, cellType2 = p.est.2)

colMeans(final_p.est)
summary(final_p.est)

hist(final_p.est$cellType1)
hist(final_p.est$cellType2)

p_0 = which(final_p.est$cellType1 < 0.01)
p_1 = which(final_p.est$cellType1 > 0.99)

length(p_0)
length(p_1)

# The End