

geneModel_X <- function (gene, d, pdDist, isoforms, lmax = length(pdDist), 
                         eLenMin = 1, verbose = 1, sam_names, mix_sams){
  
  # Status?  
  status = "OK"
  
  #--------------------------------------------------------------------------------------------------#
  # INPUT FILE CHECKING:                                                                             #
  #--------------------------------------------------------------------------------------------------#
  
  # Stops if "gene" input is not in proper list format (not a problem as call to geneModel function from
  # the isoDetector function guarantees gene is a list).  
  
  if (!is.list(gene)) {
    stop("gene must be a list\n")
  }
  if (any(!c("info", sam_names) %in% names(gene))) {
    stop("gene must be a list with components 'info' and each sample name.\n")
  }
  
  #--------------------------------------------------------------------------------------------------#
  # VARIABLE CREATION                                                                                #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # --info     : dataframe with one row per exon in the cluster detailing the "chromosome", "start", #
  #              "end", and "ExonID" for a given exon.                                               #
  # --colnms   : names of columns                                                                    #
  # --exonLens : Length of the exons in the transcript cluster.                                      #
  #--------------------------------------------------------------------------------------------------#
  
  info = gene$info
  info$exon = as.numeric(info$exon)
  colnms = c("chr", "start", "end", "exon")
  if (any(!colnms %in% names(info))) {
    stop("gene$info must be a list with components: ", colnms, 
         "\n")
  }
  exonLens = info$end - info$start
  colnms = c("count", "exons")
  
  #--------------------------------------------------------------------------------------------------#
  # INPUT FILE CHECKING                                                                              #
  #--------------------------------------------------------------------------------------------------#
  # Checks input counts, read length, fragment length distribution, and info to ensure that all      #
  # files have proper form and qualities.                                                            #
  #--------------------------------------------------------------------------------------------------#
  
  for(i in 1:length(sam_names)){
    if(any(!colnms %in% names(gene[[sam_names[i]]]))){
      stop_warn = sprintf("$%s must be a dataset with column names 'count' and 'exons'.",sam_names[i])
      stop(stop_warn)
    }
  }
  if (length(d) != 1 || (!is.numeric(d)) || d < 0) {
    stop("d must a positive number\n")
  }
  if (length(pdDist) != lmax || (!is.numeric(pdDist))) {
    stop("pdDist must be numeric vector of length lmax\n")
  }
  if (abs(sum(pdDist) - 1) > 1e-10) {
    stop("summation of pdDist must be 1\n")
  }
  if (any(pdDist < 0)) {
    stop("pdDist must be non-negative\n")
  }
  nExons = length(info$exon)
  if (!all(info$exon == 1:nExons)) {
    stop("exon IDs are not consecutive ordered numbers\n")
  }
  
  #---------------------------------------------------------------------------------------------------#
  # EXON LIST GENERATION:                                                                             #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  # exons : a list of all exon sets found in the count file.
  # wki1  : vector containing list of exon IDs used by isoform ki.
  # wki2  : vector containing all 2 combos of consecutive exons in isoform ki
  #       --e.g. Iso1 : 1,2,4,6,9,11,12,13
  #              wki2 : c("1;2","2;4","4;6","6;9","9;11","11;12","12;13")
  
  exons = as.character(NULL)
  for(i in 1:length(mix_sams)){
    exons_init = gene[[mix_sams[i]]]$exons
    exons = union(exons,exons_init)
    if (!is.null(isoforms)) {
      for (ki in 1:ncol(isoforms)) {
        isoKi = isoforms[, ki]
        wKi1 = as.character(which(isoKi == 1))
        if (length(wKi1) == 1) {
          wKi2 = character(0)
        } else {
          wKi2 = paste(wKi1[-length(wKi1)], wKi1[-1], sep = ";")
        }
        exons = union(exons,c(wKi1,wKi2))
      }
    }
  }
  # exons (at conclusion of loop):
  # Now contains all exon sets in the count file, all exons in each individual isoform, and all
  # consecutive 2 combos of exons in all known isoforms, and all exons in each sample.
  
  exons = lapply(strsplit(exons, ";"), function(v) {
    sort(as.numeric(v))
  })
  
  # exons : Now a list of same length as exons in previous iteration (all exon sets in updated)
  #         exon set list. Also, all exons are sorted in numeric order. 
  
  exn1st = sapply(exons, function(v) {
    v[1]
  })
  exn2nd = sapply(exons, function(v) {
    v[2]
  })
  exn3rd = sapply(exons, function(v) {
    v[3]
  })
  
  # Will order exon sets in numeric order up to third exon in the set (and then orders based on sum of all exons)  
  
  exnsum = sapply(exons, sum)
  exnodr = order(exn1st, exn2nd, exn3rd, exnsum, na.last = FALSE)
  exons = exons[exnodr]
  
  if (verbose > 1) {
    message(sprintf("# of isoforms=%d\n", ncol(isoforms)))
  }
  
  #--------------------------------------------------------------------------------------------------#
  # DESIGN MATRIX (X) CREATION:                                                                      #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Creates the X matrix of effective lengths for the given transcript cluster. Calls to the effLen  #
  # function established in a different program.                                                     #
  #--------------------------------------------------------------------------------------------------#
  
  # Initialize X matrix to all zeros, one row for each exon set and one column for each isoform.  
  
  X = matrix(0, nrow = length(exons), ncol = ncol(isoforms))
  rownames(X) = exons
  for (i in 1:length(exons)) {
    #-----------------------------------------------------------------------------    
    # Obtain row vector of exons in i-th exon set and determine number present
    #-----------------------------------------------------------------------------
    eSeti = exons[[i]]
    nExi = length(eSeti)
    if (all(diff(eSeti) <= 1)) {
      #-----------------------------------------------------------------------------
      # If we have consecutive exons in exon set, perform the following steps:
      #    1) compute wmat, a vector containing column indices of isoforms which
      #       use all exons in the current set
      #    2) Set i-th row of X to the same value for each of these columns.
      #-----------------------------------------------------------------------------
      wmat = which(colSums(isoforms[eSeti, , drop = FALSE]) == 
                     nExi)
      X[i, wmat] = effLen(1:nExi, exonLens[eSeti], d, pdDist, 
                          lmax)
    } else {
      #-----------------------------------------------------------------------------
      # If exons in current exon set not consecutive, perform the following steps:
      #    1) Compute eskipped (a vector containing exons which could be skipped)
      #    2) Compute sum1 (sum of all exons in current exon set)
      #    3) Compute sum0 (sum of the indicators of possible skipped exons in 
      #       all isoforms: greater than 0 implies isoform in question would have
      #       had to "skip" some exon.)
      #    4) wCase1: vector of indices of isoforms which have all of these exons
      #       but need perform no skipping.
      #    5) wCase2: vector of indices of isoforms which would have had to "skip"
      #       an exon to produce a read at the current exon set.
      #-----------------------------------------------------------------------------
      eSkipped = setdiff(min(eSeti):max(eSeti), eSeti)
      sum1 = colSums(isoforms[eSeti, , drop = FALSE])
      sum0 = colSums(isoforms[eSkipped, , drop = FALSE])
      wCase1 = which(sum1 == length(eSeti) & sum0 == 0)
      if (length(wCase1) > 0) {
        #---------------------------------------------------------------------------
        # Computes effective length for isoforms where no skipping is necessary
        #---------------------------------------------------------------------------
        X[i, wCase1] = effLen(1:nExi, exonLens[eSeti], 
                              d, pdDist, lmax)
      }
      wCase2 = which(sum1 == length(eSeti) & sum0 > 0)
      if (length(wCase2) > 0) {
        #---------------------------------------------------------------------------
        # If exon "skipping" is necessary, perform the following steps:
        #    1) Determine the unique skipping arrangements for all isoforms with
        #       regard to current exon set. (Isoforms with same exons "skipped"
        #       can be grouped together as effective length will be identical)
        #    2) w3 (loop): vector containing isoforms that match current skipping
        #       paradigm.
        #    3) w2 : index of representative isoform for current skipping paradigm.
        #---------------------------------------------------------------------------
        unqIso = unique(isoforms[eSkipped, wCase2, drop = FALSE], 
                        MARGIN = 2)
        for (jj in 1:ncol(unqIso)) {
          isojj = isoforms[eSkipped, wCase2, drop = FALSE] - 
            unqIso[, jj]
          w3 = wCase2[which(colSums(abs(isojj)) == 0)]
          w2 = w3[1]
          w2Exn = which(isoforms[, w2] > 0)
          ids = match(eSeti, w2Exn)
          wbrk = which(diff(ids) > 1)
          if (length(wbrk) <= 1) {
            rjs = exonLens[w2Exn]
            X[i, w3] = effLen(ids, rjs, d, pdDist, lmax)
          }
        }
      }
    }
  }
  
  #---------------------------------------------------------------------------------
  # w2rm Creation: Identify which columns of X (which isoforms) have effective
  #                effective lengths less than some minimal value.
  #                Used to remove these isoforms from consideration.
  #---------------------------------------------------------------------------------
  
  w2rm = which(apply(X, 2, function(v) {
    all(v <= eLenMin)
  }))
  if (length(w2rm) > 0) {
    isoforms = isoforms[, -w2rm, drop = FALSE]
    X = X[, -w2rm, drop = FALSE]
  }
  X = X + eLenMin
  if (verbose > 1) {
    message(sprintf("size of design matrix X = (%d,%d)\n", 
                    nrow(X), ncol(X)))
  }
  
  gm = list(info=info, candiIsoform = isoforms, X = X, exon_sets = exons)
  
  #--------------------------------------------------------------------------------------------------#
  # OUTPUT VECTOR CREATION - Cut                                                                     #
  #--------------------------------------------------------------------------------------------------#
  
  # creates output vector y, which is 0 if not found in count file and value find in count file if it
  # is present in the count file.
  # r = 0
  # gm = list(info=info, candiIsoform = isoforms, X = X)
  
  # for(j in 1:length(sam_names)){
  #   r = r+1
  #   count = gene[[sam_names[j]]]
  #   y = rep(0, nrow(X))
  #   exnsStr = sapply(exons, paste, collapse = ";")
  #   mm = match(count$exon, exnsStr)
  #   mn = which(mm!="NA")
  #   mm = mm[mn]
  #   y[mm] = count$count[mn]
  #   countN = data.frame(exons = exnsStr, count = y)
  #   y_out = sprintf("y_%s",sam_names[j])
  #   cn_out = sprintf("countN_%s",sam_names[j])
  #   gm[["info_status"]] = gene[["info_status"]]
  #   gm[[y_out]] = y
  #   gm[[cn_out]] = countN
  # }
  
  gm
}