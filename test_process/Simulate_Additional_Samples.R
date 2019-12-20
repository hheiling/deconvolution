# Simulate additional samples using output from Step 5

# X * alpha.est = mean of negative binomial distribution
# Need to find overdispersion paramter - use glm.nb with X*alpha.est as a covariate
# Use mean = X * alpha.est and overdispersion parameter to simulate additional count data
## Start with adding 4 more datasets (total of 5). Thoughts: 10?
# Add additional simulated data to Step 4 output

# Load the modified_sig_geneMod object from result of Step 4
step4_result = sprintf("%s/Step4_modified_sig_geneMod.RData", prefix)
load(step4_result)

# Load the pure_est object from result of Step 5
step5_result = sprintf("%s/Step5_pure_est.RData", prefix)
load(step5_result)

#-------------------------------------------------------------------#
# Experimentation                                                   #
#-------------------------------------------------------------------#
library(MASS)
samp_info = modified_sig_geneMod[[1]]$Sample_Info
cellType = samp_info$Cell_Type
cellType_pure = cellType[which(cellType != "mix")]
label = samp_info$Label
label_pure = label[which(label != "mix")]
tmp.data = pure_est[[1]]

test_data = tmp.data$chr18_1
X = test_data$X
alpha = test_data$alpha.est
sim_num = 4
i = 2
for(i in 1:length(cellType)){
  ct_name = cellType_pure[i]
  y_label = sprintf("y_%s", label_pure[i])
  y_k = test_data[[y_label]]
  X_a = X %*% alpha[,i]
  fit_pois = glm(y_k ~ X_a, family = "poisson")
  fit = glm.nb(y_k ~ X_a, link = log, control = glm.control(maxit = 50),
               etastart = (cbind(1, X_a) %*% fit_pois$coefficients))
  theta = fit$theta
  nb_means = exp((cbind(1, X_a) %*% fit$coefficients))
  y_sim = matrix(0, nrow = nrow(X), ncol = sim_num)
  for(j in 1:ncol(y_sim)){
    y_sim[,j] = rnegbin(n = nrow(X), mu = nb_means, theta = theta)
  }
  # y_sim2 = matrix(0, nrow = nrow(X), ncol = sim_num)
  # for(j in 1:ncol(y_sim)){
  #   y_sim2[,j] = rpois(n = nrow(X), lambda = nb_means)
  # }
}

#-------------------------------------------------------------------#
# Simulation                                                        #
#-------------------------------------------------------------------#
library(MASS)
library(stringr)
samp_info = modified_sig_geneMod[[1]]$Sample_Info
cellType = samp_info$Cell_Type
cellType_pure = cellType[which(cellType != "mix")]
label = samp_info$Label
label_pure = label[which(label != "mix")]
tmp.data = pure_est[[1]]

sim_num = 4

sim_data = list()

seed = 1

clusters = names(tmp.data)
for(clust in clusters){
  # cat("clust: ", clust, "\n")
  clust_data = tmp.data[[clust]]
  alpha = clust_data$alpha.est
  
  for(i in 1:length(cellType_pure)){
    # cat("i: ", i, "\n")
    ct_name = cellType_pure[i]
    y_label = sprintf("y_%s", label_pure[i])
    y = clust_data[[y_label]]
    X = clust_data$X
    X_a = X %*% alpha[,i]
    # Potential solution: fit a poisson model, then use these values as starting values
    # fit_pois = glm(y ~ X - 1, family = "poisson")
    # fit = glm.nb(y ~ X - 1, etastart = (X %*% fit_pois$coefficients), init.theta = 0.2)
    # fit_pois = glm(y ~ X_a, family = "poisson")
    # fit = glm.nb(y ~ X_a, link = log, control = glm.control(maxit = 50),
    #              etastart = (cbind(1,X_a) %*% fit_pois$coefficients),
    #              init.theta = 0.2)
    fit = glm.nb(y ~ X_a - 1, link = log)
    theta = fit$theta
    if(theta > 1){
      cat("clust: ", clust, "\n")
      cat("i: ", i, "\n")
      cat("theta: ", theta, "\n")
      theta = runif(1, 0.1, 0.3)
      cat("new theta: ", theta, "\n")
    }
    nb_means = exp(cbind(X_a) %*% fit$coefficients)
    nb_means = ifelse(nb_means > (5*max(y)), max(y), nb_means)
    y_sim = matrix(0, nrow = nrow(X), ncol = sim_num)
    set.seed(seed)
    for(j in 1:ncol(y_sim)){
      y_sim[,j] = rnegbin(n = nrow(X), mu = nb_means, theta = theta)
    }

    y_cols = str_c(sprintf("y_%s", ct_name), "_ref",2:(sim_num+1))
    colnames(y_sim) = y_cols

    seed = seed + 1

    sim_data[[clust]][[ct_name]] = list(y_sim = y_sim, seed = seed, mu = nb_means, theta = theta)
    
    if(anyNA(y_sim)){
      cat("clust: ", clust, "\n")
      cat("i: ", i, "\n")
      cat("at least one NA in y_sim \n")
    }
    
    # if(theta > 1){
    #   cat("clust: ", clust, "\n")
    #   cat("i: ", i, "\n")
    #   cat("theta: ", theta, "\n")
    # }
    
  }
}

library(stringr)
prefix = "C:/Users/hheiling/Documents/GitHub/deconvolution/test_process/test_materials/"

# save(sim_data, file = str_c(prefix, "Simulated_Counts.RData"))

# list.files(path = prefix, pattern = "_counts.txt")
count_files = c(sprintf("mf_%s_counts.txt", file), "mm9_set1_counts.txt", "mm9_set2_counts.txt")
countData = str_c(prefix, count_files)
total_cts = numeric(length(countData))
for(i in 1:length(countData)){
  countsi = read.table(countData[i], as.is = T)
  counts_col = countsi[,1]
  total_cts[i] = sum(counts_col)
}


sig_geneMod = modified_sig_geneMod[[1]]

step4_sim_data = list()

for(clust in clusters){
  clust_data = sig_geneMod[[clust]]
  sim_result = sim_data[[clust]]
  
  step4_sim_data[[clust]] = clust_data
  
  for(i in 1:length(cellType_pure)){
    ct_name = cellType_pure[i]
    label_suffix = label_pure[i]
    y = clust_data[[y_label]]
    y_sim = sim_result[[ct_name]]$y_sim
    y_cols = colnames(y_sim)
    for(j in 1:ncol(y_sim)){
      y_k = y_sim[,j]
      exons = clust_data[[sprintf("countN_%s", label_suffix)]][,1]
      step4_sim_data[[clust]][[y_cols[j]]] = y_k
      step4_sim_data[[clust]][[sprintf("countN_%s", str_remove(y_cols[j], "y_"))]] = data.frame(exons = exons, count = y_k)
      step4_sim_data[[clust]][[ct_name]][["rds_exons"]] = 
        cbind(step4_sim_data[[clust]][[ct_name]][["rds_exons"]], c(total_cts[i+1] - sum(y_k), y_k))
    }
  }
}

save(step4_sim_data, file = "test_process/test_materials/Step4_Simulated_Additions.RData")
