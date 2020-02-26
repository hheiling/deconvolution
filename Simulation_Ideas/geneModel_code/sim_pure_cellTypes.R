
# Create pure cell type reference files

library(MASS)
library(stringr)

# Load needed materials

header1 = "/home/hheiling_unc_edu/"
# header2 = "C:/Users/hheiling/Documents/GitHub/"
# Simulated output prefix
prefix_sim_out = str_c(header1,"deconvolution/Simulation_Ideas/Simulated_Output")
# prefix for simulated files
prefix_pure = str_c(header1, "Simulated_Files/Pure_Sample_Counts")
# Human materials output prefix
prefix_HM = str_c(header1,"deconvolution/Simulation_Ideas/Human_Materials")
# Sourced code location prefix
prefix_code = str_c(header1,"deconvolution/Simulation_Ideas/geneModel_code")

# Load iso_exon_1000
load(file = sprintf("%s/Iso_Exon_Params_1000.RData", prefix_sim_out))
# Load iso_exon_other
load(file = sprintf("%s/Iso_Exon_Counts_Other.RData", prefix_sim_out))
# Load geneInfo
load(file = sprintf("%s/initial_geneInfo.RData", prefix_sim_out))

# rnegbin(n = number of sample values, mu = vector of means, theta = vector of theta parameters or scalar that is recycled)

# Simulate and record exon set counts 
## count.txt files saved in folder Simulated_Output/Pure Sample Counts

n = 20

rep_pure = c(str_c("0",1:9),10:n)

files1 = str_c("CT1_ref_",rep_pure,"_counts")
files2 = str_c("CT2_ref_",rep_pure,"_counts")
files3 = str_c("CT3_ref_",rep_pure,"_counts")

# Use same seeds used in "Simulation Documentation.Rmd"
set.seed(2020)
seeds = sample(1000:9999, size = 50, replace = F)

source(sprintf("%s/sim_functions.R",prefix_code))

counts_output(exonInfo_1000 = iso_exon_1000, exonInfo_other = iso_exon_other,
              theta = geneInfo$theta,
              file_labels = c(files1, files2, files3),
              folder = prefix_pure, seed = seeds[11])



#####################################################################################################