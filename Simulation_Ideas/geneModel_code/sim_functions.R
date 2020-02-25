# Functions to simulate gene counts, isoform probability distributions, and exon set counts


#-----------------------------------------------------------------------------#
# Create gene-level output for all J cell types                                #
#-----------------------------------------------------------------------------#

# Create gene_level output for all three cell types
gene_level = function(total_cts, gene_alpha, seed){
  
  set.seed(seed)
  
  # Define variables
  n = nrow(total_cts)
  J = ncol(total_cts)
  CT_names = rep(1:J, each = n)
  ref_num = rep(1:n, times = J)
  
  gene_names = names(gene_alpha)
  
  # Matrix of counts for each gene
  count = matrix(NA, nrow = length(gene_alpha), ncol = n*J)
  
  # Sample from UCLA Dirichlet distribution
  prob = t(rdirichlet(n = n*J, alpha = gene_alpha))
  rownames(prob) = gene_names
  colnames(prob) = str_c("CT",CT_names,":","ref_",ref_num)
  
  # Using above prob vec, sample total_cts from multinomial distribution
  for(j in 1:J){
    for(i in 1:n){
      count[,i+n*(j-1)] = rmultinom(n = 1, size = total_cts[i,j], prob = prob[,i+n*(j-1)])
    }
  }
  
  rownames(count) = rownames(prob)
  colnames(count) = colnames(prob)
  
  # Vector of dispersion parameters
  # For each sample, assign a constant dispersion parameter (applies to all genes)
  ## Parameterization of interest: mu + mu^2 / theta
  theta = sample(90:120, n*J, replace = TRUE)
  names(theta) = colnames(prob)
  
  gene_output = list(p_mat = prob, ct_mat = count, gene_names = names(gene_alpha),
                     theta = theta)
  
  return(gene_output)
}



#-----------------------------------------------------------------------------#
# Select genes for differential expression and differential isoform usage     #
#-----------------------------------------------------------------------------#

diff_genes = function(CT1_counts, nTE_filtered, num_diff = 200, seed){
  
  set.seed(seed)
  
  # Specify vars
  n = ncol(CT1_counts)
  num_genes = nrow(CT1_counts)
  gene_names = rownames(CT1_counts)
  
  # Check inputs
  if(num_diff %% 4 != 0){
    # Note: will split num_diff into diff expression and diff usage. 
    ## Of diff expression genes, will split into up and down expression
    stop("Number of genes selected for differential expression must divisible by 4")
  }
  
  counts_subset = CT1_counts[which(gene_names %in% nTE_filtered$geneId),]
  genes_nT_limit = nTE_filtered$geneId[which(nTE_filtered$nT <= 11)]
  counts_subset2 = counts_subset[which(rownames(counts_subset) %in% genes_nT_limit),]
  
  # cat("Number genes with number isoforms <= 11: ", nrow(counts_subset2), "\n")
  
  # Find genes of interest from CT1 output with counts above first p20 of counts (wrt 1000 genes of interest)
  q1 = apply(counts_subset2, 2, function(x) quantile(x, probs = 0.20))
  above_q1 = ifelse(counts_subset2 > q1, 1, 0)
  
  # Select genes for differential expression
  gene_choices = rownames(counts_subset)[which(rowSums(above_q1) == n)]
  # cat("number of gene_choices after expression level and isoform number selection: ", 
  #     length(gene_choices), "\n")
  all_diff = sample(gene_choices, num_diff, replace = F)
  diffExp = sample(all_diff, num_diff/2, replace = F)
  diffUsg = all_diff[-which(all_diff %in% diffExp)]
  
  genes_diff = matrix(NA, nrow = num_genes, ncol = 2)
  genes_diff[,1] = ifelse(gene_names %in% diffExp, 1, 0)
  genes_diff[,2] = ifelse(gene_names %in% diffUsg, 1, 0)
  
  colnames(genes_diff) = c("diffExp","diffUsg")
  rownames(genes_diff) = gene_names
  
  return(genes_diff)
  
}

#-----------------------------------------------------------------------------#
# Apply fold changes to genes for differential expression for specified genes
# from diff_genes() output
#-----------------------------------------------------------------------------#

diff_exp = function(gene_counts, n, J, CT_diffExp = 2, diff_genes_mat, propUp, seed){
  
  set.seed(seed)
  
  # Specify vars
  num_genes = nrow(gene_counts)
  all_genes = rownames(gene_counts)
  
  # Checks
  if(!(CT_diffExp %in% 1:J)){
    stop("CT_diffExp must bee a number in 1 to J")
  }else if(length(CT_diffExp) > 1){
    stop("Function ony set up to handle one cell type with diff expression")
  }
  
  # Specify genes selected for differential expression
  diff_idx = which(diff_genes_mat[,"diffExp"] == 1)
  gene_choices = rownames(diff_genes_mat[diff_idx,])
  num_diffExp = sum(diff_genes_mat[,"diffExp"])
  num_upExp = round(num_diffExp*propUp)
  num_downExp = num_diffExp - num_upExp
  
  # Initialize log2 fold change matrix to be all 0
  fc = matrix(matrix(0, nrow = num_genes, ncol = 1))
  
  # Update fold change matrix
  up_diffExp = sample(gene_choices, num_upExp, replace = F)
  down_diffExp = gene_choices[-which(gene_choices %in% up_diffExp)]
  
  fc[which(all_genes %in% up_diffExp),1] = runif(n = num_upExp, min = log2(1.5), max = log2(2))
  fc[which(all_genes %in% down_diffExp),1] = runif(n = num_downExp, min = -log2(2), 
                                                   max = -log2(1.5))
  
  rownames(fc) = all_genes
  colnames(fc) = str_c("CT",CT_diffExp)
  
  # Apply fold change matrix to gene counts of cell type of interest
  
  ## Multiply CT_diffExp counts by 2^(fc + rnorm(1, mean = 0, sd = 0.05))
  gene_cts = gene_counts[,(1+n*(CT_diffExp-1)):(n*CT_diffExp)]
  fc_rand = matrix(fc, nrow = nrow(fc), ncol = n) + matrix(rnorm(n = nrow(fc)*n, mean = 0, sd = 0.05), nrow = nrow(fc))
  gene_cts_fc = gene_cts * 2^(fc_rand)
  ## Determine proportion of counts for gene_choices when no fold change
  propA = colSums(gene_cts[which(all_genes %in% gene_choices),]) / colSums(gene_cts)
  propB = colSums(gene_cts_fc[which(all_genes %in% gene_choices),]) / colSums(gene_cts_fc)
  prop_R = propA / propB
  print("proportion ratio")
  print(prop_R)
  ## Multiply ratio of proportions to the counts affected by fc
  gene_cts_fc[which(fc != 0),] = gene_cts_fc[which(fc != 0),] * 
    matrix(prop_R, nrow = num_diffExp, ncol = n, byrow = T)
  
  gene_counts_new = gene_counts
  gene_counts_new[,(1+n*(CT_diffExp-1)):(n*CT_diffExp)] = round(gene_cts_fc)
  colnames(gene_counts_new) = colnames(gene_counts)
  rownames(gene_counts_new) = rownames(gene_counts)
  
  return(gene_counts_new)
  
}

#-----------------------------------------------------------------------------#
# Determine Ending Fold Change b/w CT ref and CT j                            #
#-----------------------------------------------------------------------------#
calc_diffExp = function(gene_counts_new, gene_counts_orig, diff_genes_mat){
  
  # Define Variables
  
  # Rows associated with genes selected for differential expression
  rows = which(diff_genes_mat[,"diffExp"] == 1)
  # Original gene counts for CT j before fold change applied
  orig = gene_counts_orig[rows,]
  # New gene counts for CT j after fold change and proportion adjustment
  new = gene_counts_new[rows,]
  
  fc_indiv = new / orig
  fc_avg = rowMeans(fc_indiv)
  
  return(list(fc_indiv = fc_indiv, fc_avg = fc_avg))
  
}

#-----------------------------------------------------------------------------#
# Simulate exon set negative binomial means                                   #
#-----------------------------------------------------------------------------#

# Isoform and exon set information
# Note: gene clusters chosen so number isoforms I >= 3 for all genes
# iso_dist = "uniform": all probabilities will be close to 1/I (I = number isoforms)
#     Note: this type only uses max of alphaRange
# iso_dist = "outlier": one probability will be relatively high and the remaining prob
#     will be approx. evenly distributed among the remaining I-1 isoforms
#     Note: this type uses both min and max of alphaRange
# iso_dist = "paired": two probabilities will be relatively high and the remaining
#     probs will be approx. evenly distributed among the remaining I-2 isoforms
#     Note: this type uses both min and max of alphaRange
iso_exon_info = function(genes_info, nTE_filtered, 
                         iso_dist = rep("uniform", times = nrow(nTE_filtered)), 
                         alphaRange = c(20,50), EffLen_info, 
                         seed = seed){
  
  set.seed(seed)
  
  # names of 1,000 clusters of interest used in simulation
  clust_names = nTE_filtered$clustID
  # names of genes of interest
  gene_names = nTE_filtered$geneId
  # Number samples
  n = ncol(genes_info)
  
  # Check iso_dist
  iso_dist_options = unique(iso_dist)
  if(!all(iso_dist_options %in% c("uniform","outlier","paired"))){
    stop("iso_dist elements must be one of 'uniform', 'outlier', or 'paired'")
  }
  
  output = list()
  
  for(clust in clust_names){
    
    # name of gene associated with cluster
    gene = nTE_filtered$geneId[which(clust_names == clust)]
    # vector of counts for gene simulated in gene_level() function for n samples
    gene_ct = genes_info[which(gene_names == gene),]
    # Effective length matrix - ExI (num exon sets x num isoforms)
    X = EffLen_info[[clust]]$X
    # number isoforms
    I = ncol(X)
    # dirichlet alpha parameters for isoforms
    dir_dist = iso_dist[which(gene_names == gene)]
    if(dir_dist == "uniform"){
      alpha = rep(alphaRange[2], times = I)
    }else if(dir_dist == "outlier"){
      alpha = c(rep(alphaRange[1], times = (I-1)), alphaRange[2])
    }else if(dir_dist == "paired"){
      alpha = c(rep(alphaRange[1], times = (I-2)), rep(alphaRange[2], times = 2))
    }
    
    # isoform probability matrix - Ixn (col = isoform probability vector associated with sample i)
    rho = t(rdirichlet(n = n, alpha = alpha))
    candiIsoform = EffLen_info[[clust]]$candiIsoform
    rownames(rho) = colnames(candiIsoform)
    colnames(rho) = str_c("ref_",1:n)
    
    # scaling factor for gene
    r_g = numeric(n)
    # coefficient for mu_g = X_g %*% beta_g
    beta = matrix(NA, nrow = I, ncol = n)
    
    for(i in 1:n){
      r_g[i] = gene_ct[i] / sum(X %*% rho[,i])
      beta[,i] = rho[,i] * r_g[i]
    }
    
    # Find exon sets corresponding to rows of X
    exon_sets = EffLen_info[[clust]]$ExonSetLabels
    
    # negative binomial means for the exon sets within cluster
    # result: each col = neg bin means for sample i of n samples,
    #   each row corresponds with (possible) exon sets
    mu = X %*% beta
    rownames(mu) = exon_sets
    colnames(mu) = str_c("ref_",1:n)
    
    output[[clust]] = list(iso_alpha = alpha, rho = rho, mu = mu, exon_sets = exon_sets)
  }
  
  return(output)
  
}

#-----------------------------------------------------------------------------#
# Simulate exon set counts for 'other' genes                                  #
#-----------------------------------------------------------------------------#

# Isoform and exon set information
# Note: gene clusters chosen so number isoforms I >= 3 for all genes
# E = number of singular exon sets
# iso_dist = "uniform": all probabilities will be close to 1/E
#     Note: this type only uses max of alphaRange
# iso_dist = "outlier": one probability will be relatively high and the remaining prob
#     will be approx. evenly distributed among the remaining E-1 isoforms
#     Note: this type uses both min and max of alphaRange
other_exonset_count = function(genes_info, nTE_other, exon_sets_other, 
                               iso_dist = rep("uniform", times = nrow(nTE_other)), 
                               alphaRange = c(20,50), seed = seed){
  set.seed(seed)
  
  # names of 'other' clusters 
  clust_names = nTE_other$clustID
  # names of 'other' genes
  gene_names = nTE_other$geneId
  # Number samples
  n = ncol(genes_info)
  
  # Check iso_dist
  iso_dist_options = unique(iso_dist)
  if(!all(iso_dist_options %in% c("uniform","outlier"))){
    stop("iso_dist elements must be one of 'uniform' or 'outlier'")
  }
  
  output = list()
  
  for(clust in clust_names){
    
    # name of gene associated with cluster
    gene = gene_names[which(clust_names == clust)]
    # vector of counts for gene simulated in gene_level() function for n samples
    gene_ct = genes_info[which(gene_names == gene),]
    # singular exon sets
    exon_sets = exon_sets_other[[clust]]
    
    # Distribute gene counts to singular exon sets according to iso_dist specification
    ## E = number singular exon sets
    E = length(exon_sets)
    ## dirichlet alpha parameters
    dir_dist = iso_dist[which(gene_names == gene)]
    if(dir_dist == "uniform"){
      alpha = rep(alphaRange[2], times = E)
    }else if(dir_dist == "outlier"){
      alpha = c(rep(alphaRange[1], times = (E-1)), alphaRange[2])
    }
    
    # exon set probability matrix - Exn (col = exon set probability vector associated with sample i)
    rho = t(rdirichlet(n = n, alpha = alpha))
    rownames(rho) = exon_sets
    colnames(rho) = str_c("ref_",1:n)
    
    # Determine exon set counts from multinomial distriution
    exon_set_cts = matrix(NA, nrow = E, ncol = n)
    for(i in 1:n){
      exon_set_cts[,i] = rmultinom(n = 1, size = gene_ct[i], prob = rho[,i])
    }
    
    rownames(exon_set_cts) = exon_sets
    colnames(exon_set_cts) = colnames(genes_info)
    
    output[[clust]] = list(exon_sets = exon_sets, exon_set_cts = exon_set_cts)
  }
  
  return(output)
  
}

#-----------------------------------------------------------------------------#
# Create Pure CT Reference Count Files                                        #
#-----------------------------------------------------------------------------#

counts_output = function(exonInfo_1000, exonInfo_other, theta, file_labels,
                         folder = "Human_Materials/Pure Sample Counts", seed){
  
  set.seed(seed)
  
  # Checks
  
  
  # Define variables
  ## Number cell types
  J = length(exonInfo_1000)
  ## Number samples per cell type
  n = length(theta) / J
  
  output_1000 = list()
  
  for(ct in 1:J){
    
    ct_info = exonInfo_1000[[ct]]
    ct_theta = theta[(1 + n*(ct-1)):(n*ct)]
    
    for(clust in names(ct_info)){
      
      mu = ct_info[[clust]]$mu
      exon_sets = ct_info[[clust]]$exon_sets
      
      # Initialize ES_labels (exon set labels) and record of counts (counts_record) in first cluster
      if(names(ct_info)[1] == clust){
        
        ES_labels = exon_sets
        counts_record = matrix(NA, nrow = length(exon_sets), ncol = n)
        
        for(i in 1:n){
          counts_record[,i] = rnegbin(n = length(mu[,i]), mu = mu[,i], theta = ct_theta[i])
        }
        
      }else{ # End IF
        
        ES_labels = c(ES_labels, exon_sets)
        counts_subset = matrix(NA, nrow = length(exon_sets), ncol = n)
        
        for(i in 1:n){
          counts_subset[,i] = rnegbin(n = length(mu[,i]), mu = mu[,i], theta = ct_theta[i])
        }
        
        counts_record = rbind(counts_record, counts_subset)
        
      } # End ELSE of IF-ELSE
      
    } # End clust for-loop
    
    rownames(counts_record) = ES_labels
    
    output_1000[[ct]] = list(ES_labels = ES_labels, counts = counts_record)
    
  } # End ct for-loop
  
  output_other = list()
  
  for(ct in 1:J){
    
    ct_info = exonInfo_other[[ct]]
    
    for(clust in names(ct_info)){
      
      exon_sets = ct_info[[clust]]$exon_sets
      
      # Initialize ES_labels (exon set labels) and record of counts (counts_record) in first cluster
      if(names(ct_info)[1] == clust){
        
        ES_labels = exon_sets
        counts_record = ct_info[[clust]]$exon_set_cts
        
      }else{ # End IF
        
        ES_labels = c(ES_labels, exon_sets)
        counts_record = rbind(counts_record, ct_info[[clust]]$exon_set_cts)
        
      } # End ELSE of IF-ELSE
      
    } # End clust for-loop
    
    output_other[[ct]] = list(ES_labels = ES_labels, counts = counts_record)
    
  } # End ct for-loop
  
  for(ct in 1:J){
    
    ct_files = file_labels[(1 + n*(ct-1)):(n*ct)]
    
    counts_combo = rbind(output_1000[[ct]]$counts, output_other[[ct]]$counts)
    ES_labels_all = c(output_1000[[ct]]$ES_labels, output_other[[ct]]$ES_labels)
    
    for(i in 1:n){
      df = data.frame(counts = counts_combo[,i], exons = ES_labels_all)
      write.table(df, file = sprintf("%s/%s.txt", folder, ct_files[i]), 
                  row.names = F, col.names = F)
    } # End i for-loop
    
  } # End ct for-loop
  
}

#-----------------------------------------------------------------------------#
# Simulate Mixture Count Files                                                #
#-----------------------------------------------------------------------------#

mix_creation = function(set_mixSim, out_folder, file_labels, total_cts, probs, seed){
  
  set.seed(seed)
  
  # Define variables
  ## Number mixture replicates to create
  mix_rep = nrow(total_cts)
  ## Number cell types
  J = ncol(probs)
  ## Number pure reference samples per cell type (assume equal across all cell types)
  M = length(set_mixSim[[1]])
  
  # Checks
  if(any(rowSums(probs) != 1)){
    stop("rows of probs must add to 1")
  }
  
  # List of pure reference sample count data.frames
  df_list = list()
  for(j in 1:J){
    pure_files = set_mixSim[[j]]
    files_list = list()
    for(f in 1:length(pure_files)){
      df = read.table(file = pure_files[f], as.is = T)
      colnames(df) = c("count","exons")
      files_list[[f]] = df
    }
    df_list[[j]] = files_list
  }
  
  # exon set labels (assume same across all pure reference samples)
  exon_sets = df_list[[1]][[1]]$exons
  # Number exon sets (assume equal across all pure reference samples)
  E = length(exon_sets)
  
  for(k in 1:nrow(probs)){
    # Identify prob vector
    p = probs[k,]
    # Randomly select counts files from each cell type
    pure_counts = matrix(NA, nrow = E, ncol = J)
    for(j in 1:J){
      counts_vec = df_list[[j]][[sample(1:M, size = 1)]]$count
      pure_counts[,j] = counts_vec
    }
    
    # Calculate ratio of total counts between mixture replicate and pure reference counts
    cts_Ratio = matrix(NA, nrow = mix_rep, ncol = J)
    for(rep in 1:mix_rep){
      cts_Ratio[rep,] = total_cts[rep,k] / colSums(pure_counts)
    }
    
    # Multiply p and cts_Ratio to appropriate columns of pure_counts to get mixture sample components
    ## Round results and add results across exon sets
    mixture = list()
    for(rep in 1:mix_rep){
      mix_components = pure_counts * matrix(p, nrow = nrow(pure_counts), ncol = J, byrow = T) *
        matrix(cts_Ratio[rep,], nrow = nrow(pure_counts), ncol = J, byrow = T)
      mixture[[rep]] = rowSums(round(mix_components))
    }
    
    # Save mixture results in counts.txt files
    for(rep in 1:mix_rep){
      label = file_labels[rep,k]
      df_mix = data.frame(count = mixture[[rep]], exons = exon_sets)
      write.table(df_mix, file = sprintf("%s/%s.txt", out_folder, label), col.names = F, row.names = F)
    }
  }
  
}

mix_creation2 = function(set_mixSim, out_folder, file_labels, total_cts, probs, seed){
  
  set.seed(seed)
  
  # Define variables
  ## Number mixture replicates to create
  mix_rep = length(total_cts)
  ## Number cell types
  J = length(probs)
  ## Number pure reference samples per cell type (assume equal across all cell types)
  M = length(set_mixSim[[1]])
  
  # Checks
  if(sum(probs) != 1){
    stop("probs must add to 1")
  }
  
  # List of pure reference sample count data.frames
  df_list = list()
  for(j in 1:J){
    pure_files = set_mixSim[[j]]
    files_list = list()
    for(f in 1:length(pure_files)){
      df = read.table(file = pure_files[f], as.is = T)
      colnames(df) = c("count","exons")
      files_list[[f]] = df
    }
    df_list[[j]] = files_list
  }
  
  # exon set labels (assume same across all pure reference samples)
  exon_sets = df_list[[1]][[1]]$exons
  # Number exon sets (assume equal across all pure reference samples)
  E = length(exon_sets)
  
  # for(k in 1:nrow(probs)){
  #   # Identify prob vector
  #   p = probs[k,]
  #   # Randomly select counts files from each cell type
  #   pure_counts = matrix(NA, nrow = E, ncol = J)
  #   for(j in 1:J){
  #     counts_vec = df_list[[j]][[sample(1:M, size = 1)]]$count
  #     pure_counts[,j] = counts_vec
  #   }
  #   
  #   # Calculate ratio of total counts between mixture replicate and pure reference counts
  #   cts_Ratio = matrix(NA, nrow = mix_rep, ncol = J)
  #   for(rep in 1:mix_rep){
  #     cts_Ratio[rep,] = total_cts[rep,k] / colSums(pure_counts)
  #   }
  #   
  #   # Multiply p and cts_Ratio to appropriate columns of pure_counts to get mixture sample components
  #   ## Round results and add results across exon sets
  #   mixture = list()
  #   for(rep in 1:mix_rep){
  #     mix_components = pure_counts * matrix(p, nrow = nrow(pure_counts), ncol = J, byrow = T) *
  #       matrix(cts_Ratio[rep,], nrow = nrow(pure_counts), ncol = J, byrow = T)
  #     mixture[[rep]] = rowSums(round(mix_components))
  #   }
  #   
  #   # Save mixture results in counts.txt files
  #   for(rep in 1:mix_rep){
  #     label = file_labels[rep,k]
  #     df_mix = data.frame(count = mixture[[rep]], exons = exon_sets)
  #     write.table(df_mix, file = sprintf("%s/%s.txt", out_folder, label), col.names = F, row.names = F)
  #   }
  # }
  
  # Identify prob vector
  p = probs
  # Randomly select counts files from each cell type
  pure_counts = matrix(NA, nrow = E, ncol = J)
  for(j in 1:J){
    counts_vec = df_list[[j]][[sample(1:M, size = 1)]]$count
    pure_counts[,j] = counts_vec
  }

  # Calculate ratio of total counts between mixture replicate and pure reference counts
  cts_Ratio = matrix(NA, nrow = mix_rep, ncol = J)
  for(rep in 1:mix_rep){
    cts_Ratio[rep,] = total_cts[rep,k] / colSums(pure_counts)
  }

  # Multiply p and cts_Ratio to appropriate columns of pure_counts to get mixture sample components
  ## Round results and add results across exon sets
  mixture = list()
  for(rep in 1:mix_rep){
    mix_components = pure_counts * matrix(p, nrow = nrow(pure_counts), ncol = J, byrow = T) *
      matrix(cts_Ratio[rep,], nrow = nrow(pure_counts), ncol = J, byrow = T)
    mixture[[rep]] = rowSums(round(mix_components))
  }

  # Save mixture results in counts.txt files
  for(rep in 1:mix_rep){
    label = file_labels[rep,k]
    df_mix = data.frame(count = mixture[[rep]], exons = exon_sets)
    write.table(df_mix, file = sprintf("%s/%s.txt", out_folder, label), col.names = F, row.names = F)
  }

  
}

#-----------------------------------------------------------------------------#
# Simulate Fragment Length Distribution Files                                 #
#-----------------------------------------------------------------------------#

fragLens_out = function(total_reads, mean = 300, SD = 50, lenMin = 150, lenMax = 600,
                        out_folder, file_names, seed){
  
  set.seed(seed)
  
  # Define variables
  mix_rep = nrow(total_reads)
  num_pCombos = ncol(total_reads)
  
  for(rep in 1:mix_rep){
    
    for(p in 1:num_pCombos){
      # fragLens_dist() in geneModel_code/fragLens_dist.cpp file
      freq_dist = fragLens_dist(total_reads[rep,p], mean, SD, lenMin, lenMax)
      freq_dist = freq_dist[which(freq_dist[,1] > 0),]
      write.table(freq_dist, file = sprintf("%s/%s.txt",out_folder,file_names[rep,p]), col.names = F, row.names = F)
    }
    
  }
  
}
