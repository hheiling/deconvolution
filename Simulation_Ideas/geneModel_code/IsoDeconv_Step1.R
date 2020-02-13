

IsoDeconv_Step1 = function(directory = NULL, mix_files, pure_ref_files, fraglens_files,
                           bedFile, knownIsoforms, discrim_genes = NULL, 
                           readLen, lmax = 600, eLenMin = 1){
  
  #-----------------------------------------------------------------------------#
  # Download all needed files, convert input items to useful formats            #
  #-----------------------------------------------------------------------------#
  
  if(!is.null(directory)){
    setwd(directory)
  }
  
  countData_mix = mix_files
  countData_pure = pure_ref_files[,1] 
  
  cellTypes_pure = as.character(pure_ref_files[,2])
  # ctpure_names = pure cellTypes names
  ctpure_names = levels(as.factor(cellTypes_pure))
  
  labels_pure = character(length(countData_pure))
  
  for(type in ctpure_names){
    locations = which(cellTypes_pure == type)
    labels_type = cellTypes_pure[locations]
    num_ref = length(labels_type)
    labels_ref = str_c(labels_type, "_ref", 1:num_ref)
    labels_pure[locations] = labels_ref
  }
  
  # Download all count text files, compute total counts for each file
  # Output: list with elements total_cts, counts_list
  ## See comp_total_cts() function under "Internal isoDeconvMM Functions" heading later in this document
  pure_input = comp_total_cts(directory = directory, countData = countData_pure)
  
  mix_input = comp_total_cts(directory = directory, countData = countData_mix)
  
  fraglens_list = list()
  for(i in 1:length(fraglens_files)){
    fraglens = read.table(fraglens_files[i], as.is = T)
    if (ncol(fraglens) != 2) {
      stop(fraglens_files[i], " should have 2 columns for Freq and Len\n")
    }
    names(fraglens) = c("Freq", "Len")
    
    fraglens_list[[i]] = fraglens
  }
  
  # Package mixture counts, pure sample counts, and fragment size information together
  
  files = list()
  
  for(i in 1:length(mix_files)){
    
    mix_counts = mix_input$counts_list[[i]]
    pure_counts = pure_input$counts_list
    
    countData = pure_counts
    countData[[(length(pure_counts)+1)]] = mix_counts
    
    total_cts = c(pure_input$total_cts, mix_input$total_cts[i])
    
    fragSizeFile = fraglens_list[[i]]
    
    files[[i]] = list(countData = countData, total_cts = total_cts,
                      fragSizeFile = fragSizeFile)
    
  }
  
  # CHECKING Presence of Isoforms File:
  # IsoDeconv requires a list of known isoforms in order to model intra-sample heterogeneity. Stops program if file not present.
  # Load knownIsoforms .RData object:
  
  if(!is.null(knownIsoforms)){
    assign("isoAll", get(load(sprintf("%s/%s", directory, knownIsoforms))))
  }else{
    stop("knownIsoforms list object must be present!")
  }
  
  # Download .bed file
  bedFile_info = read.table(sprintf("%s", bedFile), sep = "\t", as.is = TRUE)
  
  bf_colNames = c("chr", "start", "end", "exon", "score", "strand")
  
  if (ncol(bedFile_info) != 6) {
    cN = paste(bf_colNames, collapse = ", ")
    stop(bedFile, " should have 6 columns: ", cN, "\n")
  }
  
  names(bedFile_info) = bf_colNames
  
  print("Finished loading supporting files")
  
  #--------------------------------------------------------------------------------#
  # Step 1
  # Calls dev_compiled_geneMod function, which is an edit of Dr. Wei Sun's
  # geneModel creation problem (after edits, now accommodates multiple cell types)
  # geneModel() altered by:                                                                                                                        
  #    Douglas Roy Wilson, Jr. 
  # Additional edits: Edited to just give info, candiIsoforms, and X for each cluster
  #--------------------------------------------------------------------------------#
  
  labels = c(labels_pure, "mix")
  cellTypes = c(cellTypes_pure, "mix")
  
  final_geneMod = list()
  
  for(j in 1:length(files)){
    
    file = files[[i]]
    
    countData = file$countData
    total_cts = file$total_cts
    fragSizeFile = file$fragSizeFile
    
    # Call dev_compiled_geneMod function
    ## See R/isoDeconv_geneModel_revised.R for dev_compiled_geneMod() code
    fin_geneMod = dev_compiled_geneMod2(countData=countData,labels = labels,total_cts = total_cts, 
                                        cellTypes=cellTypes, bedFile=bedFile_info,knownIsoforms=isoAll,
                                        fragSizeFile=fragSizeFile,readLen=readLen,lmax=lmax,
                                        eLenMin=eLenMin,discrim_genes=discrim_genes)
    
    # Perform some checks on the geneMod output
    ## See R/rem_clust.R for rem_clust() code
    # sig_geneMod = rem_clust(geneMod = fin_geneMod,co = 5,min_ind = 0)
    
    final_geneMod[[j]] = fin_geneMod
  }
  
  print("Finished creation of gene model")
  
  EffLen_info = final_geneMod[[1]]
  
  if("Sample_Info" %in% names(EffLen_info)){
    EffLen_info = EffLen_info[-which(names(EffLen_info) == "Sample_Info")]
  }
  
  #--------------------------------------------------------------------------------#
  # Attach exon set labels associated with exon sets represented by
  # the effective length design matrices
  #--------------------------------------------------------------------------------#
  
  # Exon set labels - vector in same order as X rows
  clusters = names(EffLen_info)
  
  for(clust in clusters){
    # Extract gene name of cluster
    ## Should only have one gene per cluster (due to earlier filtering of nTE object)
    gene = unique(EffLen_info[[clust]]$info$gene)
    
    if(length(gene) > 1){
      cat("cluster ", clust, "\n")
      cat("Too many genes in cluster \n")
    }
    
    # Create strings of exon set labels (format like count file output)
    sets = EffLen_info[[clust]]$exon_sets
    labels = character(length(sets))
    
    for(e in 1:length(sets)){
      e_set = sets[[e]]
      e_labs = str_c(clust, "|", gene, "|", e_set, ";")
      labels[e] = str_flatten(e_labs)
      
    }
    
    EffLen_info[[clust]][["ExonSetLabels"]] = labels
    # rownames(EffLen_info[[clust]][["X"]]) = labels
    
  }
  
  #--------------------------------------------------------------------------------#
  # Modify effective length design matrices:
  # Remove exon sets (and corresponding rows in the X matrix) that have all 0
  # If an entry has a value above 0 but less than 1, round to 1
  # Modify ExonSetLabels to reflect removal of exons (if applicable for a cluster)
  #--------------------------------------------------------------------------------#
  
  for(clust in clusters){
    
    X = EffLen_info[[clust]]$X
    exon_sets = EffLen_info[[clust]]$ExonSetLabels
    
    keep = which(rowSums(X) > 0)
    X2 = X[keep,]
    exon_sets2 = exon_sets[keep]
    
    X3 = ifelse(X2 > 0 & X2 < 1, 1, X2)
    
    EffLen_info[[clust]][["X"]] = X3
    EffLen_info[[clust]][["ExonSetLabels"]] = exon_sets2
    
  }
  
  return(EffLen_info)
  
}

