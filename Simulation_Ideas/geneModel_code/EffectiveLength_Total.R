# Cite isoform package here

#' @export
effLen <- function (ids, rjs, d, pdDist, lmax){

  # Order ids and determine number in exon set:
  ids = sort(ids)
  k = length(ids)

#-----------------------------------------------------#
# Exon Sets of Length 1                               #
#-----------------------------------------------------#
  if (k == 1) {
    eLen = effLen1(rjs[ids], d, pdDist, lmax)
    return(eLen)
  }

#-----------------------------------------------------#
# Exon sets of Length >= 2                            #
#-----------------------------------------------------#
  difs = diff(ids)
  wbrk = which(difs > 1)

  #---------------------------------------------------#
  # All Exons Consecutive                             #
  #---------------------------------------------------#
  if (length(wbrk) == 0) {
    rjs = rjs[ids]
    ids = 1:k
    #-------------------------------------------------#
    # If more than 5 exons                            #
    #-------------------------------------------------#
    # Reduce complexity of exon set by lumping together
    # smaller exons into single "Exon Conglomerates"
    if (k > 5) {
      rjAdj = rjs
      kk = 0
      jj = 1
      while (jj <= k) {
        if (rjs[jj] < 10) {
          ww = jj + max(which(cumsum(rjs[jj:k]) < 10)) -
            1
          kk = kk + 1
          rjAdj[kk] = sum(rjs[jj:ww])
          jj = ww + 1
        } else {
          kk = kk + 1
          rjAdj[kk] = rjs[jj]
          jj = jj + 1
        }
      }
      rjs = rjAdj[1:kk]
      k = kk
      ids = 1:k
    }
    eLen = effLen1(sum(rjs), d, pdDist, lmax)
    if (eLen > 0.1) {
      for (r in 1:(k - 1)) {
        setr = combinations(k, r)
        for (i in 1:nrow(setr)) {
          seti = setr[i, ]
          idsi = ids[seti]
          eLen = eLen - effLen(idsi, rjs, d, pdDist,
                               lmax)
        }
      }
    }
  } else if (length(wbrk) == 1) {
    #---------------------------------------------------#
    # Single Jump in Exon Ids                           #
    #---------------------------------------------------#
    set1 = ids[1:wbrk]
    set2 = (ids[wbrk] + 1):(ids[wbrk + 1] - 1)
    set3 = ids[-(1:wbrk)]
    gLn1 = sum(rjs[set1])
    gLn2 = sum(rjs[set2])
    gLn3 = sum(rjs[set3])
    eLen = effLen2(gLn1, gLn2, gLn3, d, pdDist, lmax)
    if (eLen > 0.1) {
      k1 = length(set1)
      k3 = length(set3)
      if (k1 == 1 && k3 == 1) {
        return(eLen)
      }
      # Reduce complexity to assume at most 3 exons
      # in either end piece.
      if (k1 >= 3) {
        rj11 = rjs[set1[1]]
        rj13 = rjs[set1[length(set1)]]
        rj12 = gLn1 - rj11 - rj13
        if (rj12 >= d - 1) {
          return(0)
        } else {
          set1N = 1:3
          rjs1N = c(rj11, rj12, rj13)
        }
      } else {
        set1N = 1:k1
        rjs1N = rjs[set1]
      }
      if (k3 >= 3) {
        rj31 = rjs[set3[1]]
        rj33 = rjs[set3[length(set3)]]
        rj32 = gLn3 - rj31 - rj33
        if (rj32 >= d - 1) {
          return(0)
        } else {
          set3N = 1:3 + 1 + length(set1N)
          rjs3N = c(rj31, rj32, rj33)
        }
      } else {
        set3N = 1:k3 + 1 + length(set1N)
        rjs3N = rjs[set3]
      }
      rjsN = c(rjs1N, gLn2, rjs3N)
      k1N = length(set1N)
      k3N = length(set3N)
      for (r in 1:k1N) {
        for (s in 1:k3N) {
          if (r == k1N && s == k3N) {
            next
          }
          setr = combinations(k1N, r = r, v = set1N)
          sets = combinations(k3N, r = s, v = set3N)
          for (i in 1:nrow(setr)) {
            for (j in 1:nrow(sets)) {
              idsij = c(setr[i, ], sets[j, ])
              eLen = eLen - effLen(idsij, rjsN, d, pdDist,
                                   lmax)
            }
          }
        }
      }
    }
  } else {
    eLen = 0
  }
  if (eLen < 0) {
    eLen = 0
  }
  if (eLen < -10) {
    stop("ah... bug... negative effective length... \n")
  }
  eLen
}

effLen1 <- function (rj, d, pdDist, lmax)
{
  eLen = 0
  if (rj >= d) {
    lRange = d:min(rj, lmax)
    eLen = sum(pdDist[lRange] * (rj + 1 - lRange))
  }
  eLen
}

effLen2 <- function (rj, rh, rk, d, pdDist, lmax)
{
  eLen = 0
  if (rj < d || rk < d || rh + 2 * d > lmax) {
    eLen = 0
  } else {
    lRange = (2 * d + rh):min(rj + rh + rk, lmax)
    ljRange = pmin(rj, lRange - rh - d) - pmax(d, lRange -
                                                 rh - rk) + 1
    eLen = sum(pdDist[lRange] * (ljRange))
  }
}


