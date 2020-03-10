
pasteUniqu = function(v){
  it = unlist(strsplit(v, ":", fixed=TRUE))
  paste(sort(unique(it)),collapse=":")
}

annoVersion = "gencode.v15"

## when writing data into text file, it may use scientific format
## when you read it into c, and using atoi. it will make mistakes
## say 97000000 is written as 9.7e+07, and c think it is 9
## options("scipen") can control write out behavior

options(scipen=20)

## combine two overlapping exons if the overlapping part
## is more than overlapPropU of either exon
overlapPU = 0.99

## split the overlapping region if the overlapping part
## is less than overlapPropL of each of the exons
overlapPL = 0.01

## otherwise create an new exon region using the overlapping part

library(data.table)
library(stringr)

# ----------------------------------------------------------------- 
# read in data
# -----------------------------------------------------------------

inf = fread(sprintf("%s.unique.exon.gtf.gz", annoVersion))
dim(inf)
inf[1:2,]

names(inf) = c("chr", "source", "feature", "start", "end", 
"score", "strand", "frame", "anno")

dim(inf)
inf[1:2,]

table(inf$chr)
table(inf$strand)

summary(inf$end - inf$start)

# ----------------------------------------------------------------- 
# check wether exons overlap
# -----------------------------------------------------------------

nn = nrow(inf)

gaps1 = inf$start[-1] - inf$end[-nn]
gaps2 = inf$start[-1] - inf$start[-nn]

pdf(sprintf("figures/%s_gaps_unique_exons.pdf", annoVersion), width=8, height=7)
par(mar=c(5,4,2,1), mfrow=c(2,1))

hist(gaps1[abs(gaps1)<200], breaks=40, ylab="Frequency", main="",
xlab="Start of the (i+1)-th exon - End of the i-th exon")
hist(gaps2[abs(gaps2)<200], breaks=40, ylab="Frequency", main="",
xlab="Start of the (i+1)-th exon - Start of the i-th exon")

dev.off()

pdf(sprintf("figures/%s_len_unique_exons.pdf", annoVersion), width=4, height=3.5)
par(mar=c(5,4,2,1))

hist(log10(inf$end - inf$start), ylab="Frequency", main="",
xlab="Exon length, log10(bp)")

dev.off()

# -------------------------------------------------------------------
# obtain gene_id, transcript_id, gene_name, and transcript_name
# -------------------------------------------------------------------

geneId = str_extract(inf$anno, '(?<=gene_id\\s")(\\S+)(?=";)')
tranId = str_extract(inf$anno, '(?<=transcript_id\\s")(\\S+)(?=";)')
geneNm = str_extract(inf$anno, '(?<=gene_name\\s")(\\S+)(?=";)')
tranNm = str_extract(inf$anno, '(?<=transcript_name\\s")(\\S+)(?=";)')

# ----------------------------------------------------------------- 
# construct infA
# ----------------------------------------------------------------- 

infA = data.frame(geneId, tranId, geneNm, tranNm, stringsAsFactors=F)
inf  = cbind(inf[,-9], infA)

dim(inf)
inf[1:2,]

# =================================================================
# identify transcript clusters
# =================================================================

# ----------------------------------------------------------------- 
# assign cluster ID, chrom by chrom
# ----------------------------------------------------------------- 

table(inf$chr)
chrs = unique(inf$chr)

geneIDs = clusterIDs = NULL

for(chr1 in chrs){
  cat(chr1, date(), "\n")
  
  wchr    = which(inf$chr==chr1)
  infChr  = inf[wchr,]
  nn      = nrow(infChr)
  od      = order(infChr$start, infChr$end)
  if(any(diff(od) < 0)){
    stop("exons are not ordered\n")
  }
  
  gIDchr  = unique(infChr$geneId)
  length(gIDchr)

  ugIds   = sort(unique(unlist(strsplit(gIDchr, split=":"))))
  length(ugIds)
  
  ## first, those genes that share an exon
  grepIt  = grep(":", gIDchr)
  gIDcom1 = gIDchr[grepIt]
  
  ## next, those genes whose exons overlap
  cond1   = (infChr$start[-1] - infChr$end[-nn]) <= 0
  cond2   = infChr$geneId[-1] != infChr$geneId[-nn]
  
  wOverlp = which(cond1 & cond2)
  
  gIDcom2 = NULL
  
  if(length(wOverlp) > 0){
    gIDcom2 = rep("", length(wOverlp))
    for(j in 1:length(wOverlp)){
      wj = wOverlp[j]
      gIDcom2[j] = pasteUniqu(infChr$geneId[wj:(wj+1)])
    }
  }
  
  ## now, combine them
  gIDcom = unique(c(gIDcom1, gIDcom2))
  
  if(length(gIDcom) == 0){
    cut1 = 1:length(ugIds)
  }else{
    ## identify clusters by hclust using single link function
    connect = matrix(0, nrow=length(ugIds), ncol=length(ugIds))

    for(g1 in gIDcom){
      g2  = unlist(strsplit(g1, split=":"))
      wg2 = match(g2, ugIds)
      c2s  = combn(wg2, 2)
      for(k in 1:ncol(c2s)){
        c2i = c2s[1,k]
        c2j = c2s[2,k]
        connect[c2i,c2j] = 1
        connect[c2j,c2i] = 1
      }
    }
    
    diag(connect) = 1
    dM   = as.dist(1 - connect)
    h1   = hclust(dM, method="single")
    cut1 = cutree(h1, h=0.5)
    
    rm(connect, dM, h1)
  }
  
  geneIDs    = c(geneIDs, ugIds)
  clusterIDs = c(clusterIDs, paste(chr1, cut1, sep="_"))
  
}

geneId1 = sapply(strsplit(inf$geneId, split=":"), function(v){v[1]})
mat1    = match(geneId1, geneIDs)

if(any(is.na(mat1))){
  stop("something is NA here\n")  
}

all(geneIDs[mat1] == geneId1)
inf$clustId = clusterIDs[mat1]

# ----------------------------------------------------------------- 
# cluster by cluster, identify non-overlap exons
# ----------------------------------------------------------------- 

# assume the number of non-overlapping exons is no larger than 
# 5 times of the number of unique exons
infN  = matrix("", nrow=nrow(inf)*5, ncol=14)

ucIds = unique(inf$clustId)
length(ucIds)

kk   = 0
idx1 = idx2 = 1

for(uId1 in ucIds){
  
  kk = kk + 1
  if(kk %% 1000 == 0){
    cat(".")
    if(kk %% 10000 == 0){
      cat(".\n")
      cat(kk, date(), "\n")
    }
  }
  
  ## inf1 is the annotation for one cluster ID
  inf1 = inf[which(inf$clustId == uId1), ,drop=FALSE]  
  nn   = nrow(inf1)
  
  ## if this cluster only has one exon
  if(nn == 1) { 
    inf1$exonId = 1
    infN[idx1,] = as.matrix(inf1)
    idx1 = idx1 + 1
    next 
  }
  
  ## calcualte gaps
  gaps1   = inf1$start[-1] - inf1$end[-nn]  
  w2check = which(gaps1 <= 0)
    
  ## if there is no overlapping exons
  if(length(w2check) == 0){
    inf1$exonId = 1:nn
    idx2 = idx1 + nrow(inf1) - 1
    infN[idx1:idx2,] = as.matrix(inf1)
    idx1 = idx2 + 1
    next 
  }
    
  while(length(w2check) > 0){
    
    pos1 = w2check[1]
    pos2 = pos1+1
    
    ## well, we start with unique exons, but after a while...
    if(inf1$end[pos1] == inf1$end[pos2] && inf1$start[pos1] == inf1$start[pos2]){
      ## if two exons share both the start and the end poistions
      inf1$tranId[pos1] = pasteUniqu(inf1$tranId[pos1:pos2])
      inf1$tranNm[pos1] = pasteUniqu(inf1$tranNm[pos1:pos2])
      inf1$geneId[pos1] = pasteUniqu(inf1$geneId[pos1:pos2])
      inf1$geneNm[pos1] = pasteUniqu(inf1$geneNm[pos1:pos2])
      inf1              = inf1[-pos2,]
      
    }else if(inf1$end[pos1] == inf1$end[pos2]){
      ## if two exons share the end poistions

      ## if the start position difference is no greater than 2, 
      ## combine them
      if(abs(inf1$start[pos1] - inf1$start[pos2]) < 3){
        inf1$tranId[pos1] = pasteUniqu(inf1$tranId[pos1:pos2])
        inf1$tranNm[pos1] = pasteUniqu(inf1$tranNm[pos1:pos2])
        inf1$geneId[pos1] = pasteUniqu(inf1$geneId[pos1:pos2])
        inf1$geneNm[pos1] = pasteUniqu(inf1$geneNm[pos1:pos2])
        inf1              = inf1[-pos2,]
      }else{
        inf1$end[pos1]    = inf1$start[pos2] - 1
        inf1$tranId[pos2] = pasteUniqu(inf1$tranId[pos1:pos2])
        inf1$tranNm[pos2] = pasteUniqu(inf1$tranNm[pos1:pos2])
        inf1$geneId[pos2] = pasteUniqu(inf1$geneId[pos1:pos2])
        inf1$geneNm[pos2] = pasteUniqu(inf1$geneNm[pos1:pos2])
      }
    }else if(inf1$start[pos1] == inf1$start[pos2]){
      ## if two exons share the start poistions
      
      ## if the end position difference is no greater than 2, 
      ## combine them
      if(abs(inf1$end[pos1] - inf1$end[pos2]) < 3){
        inf1$tranId[pos1] = pasteUniqu(inf1$tranId[pos1:pos2])
        inf1$tranNm[pos1] = pasteUniqu(inf1$tranNm[pos1:pos2])
        inf1$geneId[pos1] = pasteUniqu(inf1$geneId[pos1:pos2])
        inf1$geneNm[pos1] = pasteUniqu(inf1$geneNm[pos1:pos2])
        inf1              = inf1[-pos2,]
      }else{
        inf1$start[pos2]  = inf1$end[pos1] + 1
        inf1$tranId[pos1] = pasteUniqu(inf1$tranId[pos1:pos2])
        inf1$tranNm[pos1] = pasteUniqu(inf1$tranNm[pos1:pos2])
        inf1$geneId[pos1] = pasteUniqu(inf1$geneId[pos1:pos2])
        inf1$geneNm[pos1] = pasteUniqu(inf1$geneNm[pos1:pos2])
      }
    }else if(inf1$end[pos2] < inf1$end[pos1]){
      ## if the 2nd exon is within the first exon
      newExon           = inf1[pos1,,drop=FALSE]
      newExon$start     = inf1$end[pos2]+1
      newExon$end       = inf1$end[pos1]
                
      inf1$end[pos1]    = inf1$start[pos2] - 1
      
      inf1$tranId[pos2] = pasteUniqu(inf1$tranId[pos1:pos2])
      inf1$tranNm[pos2] = pasteUniqu(inf1$tranNm[pos1:pos2])
      inf1$geneId[pos2] = pasteUniqu(inf1$geneId[pos1:pos2])
      inf1$geneNm[pos2] = pasteUniqu(inf1$geneNm[pos1:pos2])

      nn = nrow(inf1)
      ## new exon should be after pos2, but inf1 will be re-ordered anyway
      inf1 = rbind(inf1[1:pos1,,drop=FALSE], newExon, inf1[pos2:nn,,drop=FALSE])
      
    }else{
      
      newExon          = inf1[pos1,,drop=FALSE]
      newExon$start    = inf1$start[pos2]
      newExon$end      = inf1$end[pos1]
      newExon$tranId   = pasteUniqu(inf1$tranId[pos1:pos2])
      newExon$tranNm   = pasteUniqu(inf1$tranNm[pos1:pos2])
      newExon$geneId   = pasteUniqu(inf1$geneId[pos1:pos2])
      newExon$geneNm   = pasteUniqu(inf1$geneNm[pos1:pos2])

      inf1$end[pos1]   = newExon$start - 1
      inf1$start[pos2] = newExon$end + 1
      
      nn = nrow(inf1)
      inf1 = rbind(inf1[1:pos1,,drop=FALSE], newExon, inf1[pos2:nn,,drop=FALSE])
      
    }
    
    nn = nrow(inf1)
    if(nn==1){ break } 
    
    if(any(inf1$end < inf1$start)){
      stop("something is wrong...\n")
    }
    
    inf1 = inf1[order(inf1$start, inf1$end),]
    
    gaps1 = inf1$start[-1] - inf1$end[-nn]  
    w2check = which(gaps1 <= 0)
  }
  
  inf1$exonId = 1:nn

  idx2 = idx1 + nrow(inf1) - 1
  infN[idx1:idx2,] = as.matrix(inf1)
  idx1 = idx2 + 1

}

infN[(idx1-1):idx1,]
infN = infN[1:(idx1-1),]

infN = as.data.frame(infN, stringsAsFactors=FALSE)
names(infN) = c(names(inf), "exonId")

infN$exonId = as.numeric(infN$exonId)

dim(infN)
infN[1:2,]

# ----------------------------------------------------------------- 
# construct the anno column
# ----------------------------------------------------------------- 

geneId2use  = paste("gene_id \"", infN$geneId, "\";", sep="")
tranId2use  = paste("transcript_id \"", infN$tranId, "\";", sep="")
geneNm2use  = paste("gene_name \"", infN$geneNm, "\";", sep="")
tranNm2use  = paste("transcript_name \"", infN$tranNm, "\";", sep="")
exonId2use  = paste("exon_id \"", infN$exonId, "\";", sep="")
clustId2use = paste("clustId \"", infN$clustId, "\";", sep="")

infN$anno   = paste(geneId2use, tranId2use, geneNm2use, 
                    tranNm2use, exonId2use, clustId2use, sep=" ")

dim(infN)
infN[1:2,]

# ----------------------------------------------------------------- 
# double check the gaps between exons
# ----------------------------------------------------------------- 

nn = nrow(infN)
infN$start = as.numeric(infN$start)
infN$end   = as.numeric(infN$end)

gaps1 = infN$start[-1] - infN$end[-nn]
gaps2 = infN$start[-1] - infN$start[-nn]

w2ck = which(infN$clustId[-1] == infN$clustId[-nn])

gaps1 = gaps1[w2ck]
gaps2 = gaps2[w2ck]

summary(gaps1)
summary(gaps2)


pdf(sprintf("figures/%s_gaps_nonoverlap_exons.pdf", annoVersion), width=8, height=7)
par(mar=c(5,4,2,1), mfrow=c(2,1))

hist(gaps1[abs(gaps1)<200], breaks=40, ylab="Frequency", main="",
     xlab="Start of the (i+1)-th exon - End of the i-th exon")
hist(gaps2[abs(gaps2)<200], breaks=40, ylab="Frequency", main="",
     xlab="Start of the (i+1)-th exon - Start of the i-th exon")

dev.off()

len = infN$end - infN$start
summary(len)
table(len==0)

pdf(sprintf("figures/%s_len_nonoverlap_exons.pdf", annoVersion), width=4, height=3.5)
par(mar=c(5,4,2,1))

hist(log10(len), ylab="Frequency", main="", xlab="Exon length, log10(bp)")

dev.off()

# -----------------------------------------------------------------
# write out gtf file
# -----------------------------------------------------------------

dim(infN)
infN[1:2,]

outFile = sprintf("%s.nonoverlap.exon.gtf", annoVersion)

write.table(infN[, c(1:8,ncol(infN))], col.names = FALSE, append = FALSE, 
  file = outFile, quote = FALSE, sep = "\t", row.names = FALSE)

system(paste("gzip -f", outFile))

# -----------------------------------------------------------------
# write out bed file
# -----------------------------------------------------------------

names = paste(infN$clustId, infN$geneId, infN$exonId, sep="|")
score = rep("666", length(names))
bed   = cbind(infN$chr, infN$start, infN$end, names, score, infN$strand)

dim(bed)
bed[1:2,]

outFile = sprintf("%s.nonoverlap.exon.bed", annoVersion)

write.table(bed, col.names = FALSE, append = FALSE, 
            file = outFile, quote = FALSE, sep = "\t", row.names = FALSE)

sessionInfo()
gc()
q(save="no")


