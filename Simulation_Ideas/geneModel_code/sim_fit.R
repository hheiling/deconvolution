# Fit IsoDeconvMM model fit to simulated data

# Load IsoDeconvMM library:
# library(remotes)
# install_github("hheiling/IsoDeconvMM", force = TRUE)
library(IsoDeconvMM)
library(stringr)

# Number cell types
J = 3

header1 = "/home/hheiling_unc_edu/"
header2 = "C:/Users/hheiling/Documents/GitHub/"
# Simulated output prefix
prefix_sim_out = str_c(header1,"deconvolution/Simulation_Ideas/Simulated_Output")
# Human materials output prefix
prefix_HM = str_c(header1,"deconvolution/Simulation_Ideas/Human_Materials")
# Sourced code location prefix
prefix_code = str_c(header1,"deconvolution/Simulation_Ideas/geneModel_code")
# Simulated pure cell type files prefix
prefix_pure = str_c(header1,"Simulated_Files/Pure_Sample_Counts")
# Simulated mixture cell type files prefix
prefix_mix = str_c(header1,"Simulated_Files/Mixture_Sample_Counts")
# prefix_mix = str_c(header2,"deconvolution/Simulation_Ideas/Simulated_Output/Mixture Sample Counts")

# Source mixture sample creation function
source(sprintf("%s/sim_functions.R", prefix_code))

# Find 1,000 'genes of interest'
## load nTE_filtered data.frame
load(sprintf("%s/Homo_sapiens.GRCh37.66.nTE.filtered.RData", prefix_HM))
genes_1000 = nTE_filtered$geneId
# Find 100 genes with differential expression
## load genes_df data.frame
load(sprintf("%s/genes_simulated_w_diffExp_diffUsg.RData", prefix_sim_out))
genesDiffExp = genes_df$diffExp
# Find 100 genes with differential isoform usage
genesDiffUsage = genes_df$diffUsg

# Find pure reference samples
pure_names_full = list.files(path = prefix_pure,
                             pattern = "counts.txt$", full.names = TRUE)

cat("Number pure files found: ", length(pure_names_full), "\n")

n = length(pure_names_full) / J

CT1_files = pure_names_full[str_detect(pure_names_full, "CT1_")]
CT2_files = pure_names_full[str_detect(pure_names_full, "CT2_")]
CT3_files = pure_names_full[str_detect(pure_names_full, "CT3_")]

# Use replicates 1 through 10 for IsoDeconvMM fit algorithm
CT1_files_ref = CT1_files[1:(n/2)]
CT2_files_ref = CT2_files[1:(n/2)]
CT3_files_ref = CT3_files[1:(n/2)]

CT_ref_files = c(CT1_files_ref, CT2_files_ref, CT3_files_ref)

# Check correct selection of files 01 - 10
print(CT_ref_files)

pure_ref_files = cbind(CT_ref_files, c(rep("CT1",times=n/2),rep("CT2",times=n/2),rep("CT3",times=n/2)))

cat("dim pure_ref_files: ",dim(pure_ref_files),"\n")

# Use replicates 11 through 20 for mixture file creation
CT1_files_mixSim = CT1_files[(n/2 + 1):n]
CT2_files_mixSim = CT2_files[(n/2 + 1):n]
CT3_files_mixSim = CT3_files[(n/2 + 1):n]

set_mixSim = list()
set_mixSim[["CT1"]] = CT1_files_mixSim
set_mixSim[["CT2"]] = CT2_files_mixSim
set_mixSim[["CT3"]] = CT3_files_mixSim

# Specify file for fragment length distribution
## For this simulation, use one fragment length distribution file for each mixture file
fraglens_file = list.files(path = prefix_sim_out, pattern = "fraglens", full.names = T)

# Load probability combinations: p_combos
load(file = sprintf("%s/Probability_Combinations_Sim.RData", prefix_sim_out))

cat("Num prob combos: ",nrow(p_combos), "\n")

colnames(p_combos) = c("CT1","CT2","CT3")
rownames(p_combos) = str_c("pc_", c(str_c("0",1:9),10:nrow(p_combos)))

pc_labels = rownames(p_combos)

# set seeds
set.seed(1269)
seeds = sample(1000:9999, size = 3, replace = F)


# Specify total read counts for mixture replicates
set.seed(seeds[3])
total_cts_mix = rnorm(n = nrow(p_combos), mean = 7*10^6, sd = 10^6)
names(total_cts_mix) = pc_labels

# Create mixture samples
mix_labels = str_c("Mix_ProbCombo_",str_remove(pc_labels,"pc_"),"_counts")

cat("length of mix_labels:",length(mix_labels),"\n")
cat("First mix_label:",mix_labels[1],"\n")

batch_cut = list(1:18,19:nrow(p_combos))

# Run simulations

# Simulation 1: No initPts specified (default initPts = 1/J for all cell types)

SimResults = list()
SimSummary = list()
EffLen_Mat = list()
SimFull = list()

cellTypes = c("CT1","CT2","CT3")

for(batch in 1:length(batch_cut)){
  
  pc_nums = batch_cut[[batch]]
  
  # Probability combination label
  pc = pc_labels[pc_nums]
  
  # Create mixture files
  mix_creation2(set_mixSim = set_mixSim, out_folder = prefix_mix,
               file_labels = mix_labels[pc_nums], total_cts = total_cts_mix[pc_nums], 
               probs = p_combos[pc_nums,], seed = seeds[batch])
  
  mix_files = list.files(path = prefix_mix, pattern = "counts.txt",
                         full.names = T)
  
  mix_names = str_remove(mix_labels[pc_nums], "_counts")
  
  # Fit IsoDeconvMM algorithm
  SimResults_Full = IsoDeconvMM(directory = NULL, mix_files = mix_files,
                                pure_ref_files = pure_ref_files,
                                fraglens_files = fraglens_file,
                                bedFile = sprintf("%s/Homo_sapiens.GRCh37.66.nonoverlap.exon.bed", prefix_HM),
                                knownIsoforms = sprintf("%s/Homo_sapiens.GRCh37.66.nonoverlap.exon.knownIsoforms.RData", prefix_HM),
                                discrim_genes = union(genesDiffExp, genesDiffUsage),
                                readLen = 75, lmax = 600, eLenMin = 1,
                                mix_names = mix_names, initPts = NULL,
                                optim_options = optimControl(simple.Init = FALSE))
  
  # Save all results for reference
  SimFull[[batch]] = SimResults_Full
  
  # Save subset of results
  for(m in mix_names){
    
    results = SimResults_tmp[[m]]
    clust_names = names(results)
    
    for(clust in clust_names){
      
      SimResults[[m]][[clust]][["mix"]][["gamma.est"]] = results[[clust]][["mix"]][["gamma.est"]]
      SimResults[[m]][[clust]][["mix"]][["tau.est"]] = results[[clust]][["mix"]][["tau.est"]]
      SimResults[[m]][[clust]][["mix"]][["p.est"]] = results[[clust]][["mix"]][["p.est"]]
      SimResults[[m]][[clust]][["l_tilde"]] = results[[clust]][["l_tilde"]]
      SimResults[[m]][[clust]][["alpha.est"]] = results[[clust]][["alpha.est"]]
      SimResults[[m]][[clust]][["beta.est"]] = results[[clust]][["beta.est"]]
      SimResults[[m]][[clust]][[cellTypes[1]]][["tau.hat"]] = results[[clust]][[cellTypes[1]]][["tau.hat"]]
      SimResults[[m]][[clust]][[cellTypes[2]]][["tau.hat"]] = results[[clust]][[cellTypes[2]]][["tau.hat"]]
      SimResults[[m]][[clust]][[cellTypes[3]]][["tau.hat"]] = results[[clust]][[cellTypes[3]]][["tau.hat"]]
      SimResults[[m]][[clust]][[cellTypes[1]]][["gamma.hat"]] = results[[clust]][[cellTypes[1]]][["gamma.hat"]]
      SimResults[[m]][[clust]][[cellTypes[2]]][["gamma.hat"]] = results[[clust]][[cellTypes[2]]][["gamma.hat"]]
      SimResults[[m]][[clust]][[cellTypes[3]]][["gamma.hat"]] = results[[clust]][[cellTypes[3]]][["gamma.hat"]]
      SimResults[[m]][[clust]][["CellType_Order"]] = results[[clust]][["CellType_Order"]]
      SimResults[[m]][[clust]][["WARN"]] = results[[clust]][["WARN"]]
      
      EffLen_Mat[[m]][[clust]][["X"]] = results[[clust]][["X"]]
      EffLen_Mat[[m]][[clust]][["X.prime"]] = results[[clust]][["X.prime"]]
      
    }
    
  }
  
  # Find probability estimates for mixtures
  SimSummary[[batch]] = Summarize_Report(SimResults_Full, plots_options = plotsControl(plots = FALSE))
  
  # Delete mixture files to save space on storage
  system(sprintf("rm %s/*_counts.txt", prefix_mix))
  
}

save(SimResults, file = sprintf("%s/IsoDeconvMM_Sim_Results_200GenesP100.RData", prefix_sim_out))
save(SimSummary, file = sprintf("%s/IsoDeconvMM_Sim_Summary_200GenesP100.RData", prefix_sim_out))
save(EffLen_Mat, file = sprintf("%s/EffLen_Mat_Record_Sim_200GenesP100.RData", prefix_sim_out))

##########################################################################################################

