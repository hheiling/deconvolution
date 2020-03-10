
## when writing data into text file, it may use scientific format
## when you read it into c, and using atoi. it will make mistakes
## say 97000000 is written as 9.7e+07, and c think it is 9
## options("scipen") can control write out behavior

options(scipen=20)
annoVersion = "gencode.v15"

library(data.table)
library(stringr)

# ---------------------------------------------------------------------
# read gtf file
# ---------------------------------------------------------------------

date()
inf = fread("gencode.v15.annotation.gtf.gz")
date()

dim(inf)
inf[1:2,]

names(inf) = c("chr", "source", "feature", "start", "end", 
               "score", "strand", "frame", "anno")

sapply(c("chr", "source", "feature", "score", "strand", "frame"),
       function(xx) table(inf[,..xx], useNA="ifany"))

inf = inf[which(inf$feature == "exon"),]
dim(inf)
inf[1:2,]
inf[which(inf$source=="ENSEMBL")[1:2],]

gc()


# ---------------------------------------------------------------------
# obtain gene_id, transcript_id, gene_name, and transcript_name
# ---------------------------------------------------------------------

geneId = str_extract(inf$anno, '(?<=gene_id\\s")(\\S+)(?=";)')
tranId = str_extract(inf$anno, '(?<=transcript_id\\s")(\\S+)(?=";)')
geneNm = str_extract(inf$anno, '(?<=gene_name\\s")(\\S+)(?=";)')
tranNm = str_extract(inf$anno, '(?<=transcript_name\\s")(\\S+)(?=";)')

exonId = paste(inf$chr, inf$start, inf$end, sep=":")

tid = table(exonId)
table(tid)
sort(tid, decreasing=TRUE)[1:3]

pasteUniqu = function(v){paste(unique(v),collapse=":")}
geneId2use = tapply(geneId, exonId, pasteUniqu)

message("there are ", length(geneId2use), " unique exons.")

xx = grep(":", geneId2use)
if(length(xx) > 0){
  message(length(xx), " exons belong to more than one gene.")
  geClusters = strsplit(geneId2use[xx], split=":")
  t1         = table(sapply(geClusters, length))
  message("their distributuion is")
  print(t1)
}

tranId2use = tapply(tranId, exonId, pasteUniqu)
xx = grep(":", tranId2use)
if(length(xx) > 0){
  message(length(xx), " exons belong to more than one transcript.")
  trClusters = strsplit(tranId2use[xx], split=":")
  t1         = table(sapply(trClusters, length))
  message("their distributuion is")
  print(t1)
}

geneNm2use = tapply(geneNm, exonId, pasteUniqu)
tranNm2use = tapply(tranNm, exonId, pasteUniqu)

# ---------------------------------------------------------
# drop duplicated exons
# ---------------------------------------------------------

infNew = list()
nms = names(inf)
nms[1:8]

for(i in 1:8){
  nm1 = nms[i]  
  cat(i, nm1, "\n")

  if(nm1 == "source" || nm1 == "strand"){
    it1 = tapply(inf[[nm1]], exonId, pasteUniqu)
  }else{
    it1 = tapply(inf[[nm1]], exonId, unique)
  }
  
  if(mode(it1) == "list") { stop("hm... non unique ", nm1, "\n") }
  
  infNew[[nm1]] = it1
}

geneId2use = paste("gene_id \"", geneId2use, "\";", sep="")
tranId2use = paste("transcript_id \"", tranId2use, "\";", sep="")
geneNm2use = paste("gene_name \"", geneNm2use, "\";", sep="")
tranNm2use = paste("transcript_name \"", tranNm2use, "\";", sep="")

infNew$anno = paste(geneId2use, tranId2use, geneNm2use, tranNm2use, sep=" ")

infNew = as.data.frame(infNew)
dim(infNew)
infNew[1:2,]

# --------------------------------------------------------- 
# sort the exons
# ---------------------------------------------------------

od     = order(infNew$chr, infNew$start, infNew$end)
any(diff(od) < 0)
infNew = infNew[od,]

# ---------------------------------------------------------
# write out
# ---------------------------------------------------------

id  = paste(infNew$chr, infNew$start, infNew$end, sep=":")
uid = unique(id)

length(id)
length(uid)

outFile = sprintf("%s.unique.exon.gtf", annoVersion)

write.table(infNew, file = outFile, append = FALSE, 
  quote = FALSE, sep = "\t", row.names = FALSE, 
  col.names = FALSE)

system(paste("gzip -f", outFile))

sessionInfo()
gc()
q(save="no")
