
library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(stringr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
args

if (length(args)==0) {
  message("no argument is provided, using defaults\n")
  dataset  = "EGAD00001002671"
  sam_name = "EGAF00001331297"
  gene_anno_dir  = "_prepare_gene_anno"
  gene_anno_file = "exon_by_genes_gencode.v15.GRCh37.rds"
  bed_file       = "gencode.v15.nonoverlap.exon.bed"
} else if(length(args)==2) {
  message("two argument are provided, assume they are dataset and sam_name\n")
  # dataset  = args[1]
  # sam_name = args[2]
  eval(parse(text=args[1]))
  eval(parse(text=args[2]))
  gene_anno_dir  = "_prepare_gene_anno"
  gene_anno_file = "exon_by_genes_gencode.v15.GRCh37.rds"
  bed_file       = "gencode.v15.nonoverlap.exon.bed"
}else if(length(args)==5){
  # dataset  = args[1]
  # sam_name = args[2]
  # gene_anno_dir  = args[3]
  # gene_anno_file = args[4]
  # bed_file       = args[5]
  for(k in 1:5){
    eval(parse(text=args[k]))
  }
}else{
  stop("unexpected number of arguments")
}

# workDir   = "/Users/wsun/research/data/EGA"
# resultDir = file.path(workDir, paste0(dataset, "_result"))

workDir   = "/fh/scratch/delete90/sun_w/plittle/CS_eQTL/s5_EGA"
resultDir = file.path("/fh/scratch/delete90/sun_w/EGA", paste0(dataset, "_result"))

readLen   = 100

gene_anno_file = file.path(gene_anno_dir, gene_anno_file)
bed_file = file.path(gene_anno_dir, bed_file)

print(sprintf("sam_name: %s", sam_name))

# ------------------------------------------------------------------------
# read in sample information
# ------------------------------------------------------------------------

meta = readRDS("data/blueprint_meta_info.rds")
dim(meta)
meta[1:2,]

table(meta$LIBRARY_LAYOUT, meta$DATASET_ID)

w2do = which(meta$EGAF == sam_name)
w2do

bam_files = list.files(file.path(workDir, dataset, sam_name), pattern=".bam$")
bam_files
meta$filename[w2do]

mat1 = str_detect(meta$filename[w2do], bam_files)
if(sum(mat1) != 1){ stop("non-unique match") }

bam_file = file.path(workDir, dataset, sam_name, bam_files[mat1])
bam_filtered   = gsub(".bam$", "_filtered.bam", bam_file)

bam_file
singleEnd = meta$LIBRARY_LAYOUT[w2do]
singleEnd

# ------------------------------------------------------------------------
# counting
# ------------------------------------------------------------------------

ct1 = countBam(bam_file)
print("done with first counting!\n")

ct1$nucleotides/ct1$records

if(abs(ct1$nucleotides/ct1$records - readLen) > 5){
  stop("looks like readLen is not expected!")
}

# ------------------------------------------------------------------------
# index bam file if needed
# ------------------------------------------------------------------------

if(! file.exists(paste0(bam_file, ".bai"))){
  indexBam(bam_file)
}

# ------------------------------------------------------------------------
# getUnique and filtering
# These RNA-seq reads were mapped by STAR, and for STAR, The mapping 
# quality MAPQ (column 5) is 255 for uniquely mapping reads, 
# and int(-10*log10(1-1/Nmap)) for multi-mapping reads
# ------------------------------------------------------------------------

if(singleEnd=="SINGLE"){properPair = FALSE} else{properPair = TRUE}

flag1  = scanBamFlag(isUnmappedQuery=FALSE, isSecondaryAlignment=FALSE,
                     isDuplicate=FALSE, isNotPassingQualityControls=FALSE,
                     isSupplementaryAlignment=FALSE, isProperPair=properPair)

param1 = ScanBamParam(flag=flag1, what="seq", mapqFilter=255)

filterBam(bam_file, destination=bam_filtered, param=param1)
print("done with filtering!")

# ------------------------------------------------------------------------
# counting again
# ------------------------------------------------------------------------

ct2 = countBam(bam_filtered)
print("done with second counting!\n")

print("the total number of reads/nucleotides before/after filtering:")
print(ct1)
print(ct2)

ct2$nucleotides/ct2$records

# ------------------------------------------------------------------------
# calculate total read count (TReC) per gene
# ------------------------------------------------------------------------

if(singleEnd=="SINGLE"){singleEnd_i = TRUE} else{singleEnd_i = FALSE}

genes   = readRDS(gene_anno_file)
bamfile = BamFileList(bam_filtered, yieldSize=1000000)

se = summarizeOverlaps(features=genes, reads=bamfile, mode="Union",
                       singleEnd=singleEnd_i, ignore.strand=TRUE)

ct = as.data.frame(assay(se))

print("done with TReC!")

write.table(ct, file = file.path(resultDir, sprintf("%s_trec.txt", sam_name)), 
            quote = FALSE, sep = "\t", eol = "\n")

# ------------------------------------------------------------------------
# extract reads per exon-set
# ------------------------------------------------------------------------

un_exons = fread(bed_file)
dim(un_exons)
un_exons[1:2,]
names(un_exons) = c("chr", "start", "end", "info", "score", "strand")
un_exons = data.frame(un_exons, stringsAsFactors = FALSE)

table(un_exons$strand)
un_exons = un_exons[which(un_exons$strand %in% c("+", "-")),]
rownames(un_exons) = un_exons$info
dim(un_exons)
un_exons[1:2,]

un_exons_GR = makeGRangesFromDataFrame(un_exons)

se_exon = summarizeOverlaps(features=un_exons_GR, reads=bamfile, mode="Union",
                       singleEnd=singleEnd_i, ignore.strand=TRUE)

ct_exon = as.data.frame(assay(se_exon))
dim(ct_exon)
ct_exon[1:2,,drop=FALSE]

print("done with TReC per exon-set!")

write.table(ct_exon, file = file.path(resultDir, sprintf("%s_trec_exon_set.txt", sam_name)), 
            quote = FALSE, sep = "\t", eol = "\n")

# ------------------------------------------------------------------------
# summarize fragment length distribution
# ------------------------------------------------------------------------

if(!singleEnd_i){
  tempFile    = file.path(workDir, dataset, sam_name, paste0(sam_name, "_temp.txt"))
  fragLenFile = file.path(resultDir, paste0(sam_name, "_fragLen.txt"))
  
  cmd1 = sprintf("samtools view -f 65 %s | awk '{print ($8>=$4) ? $8-$4+%s : $4-$8+%s}' > %s",
                 bam_filtered, readLen, readLen, tempFile)
  system(cmd1)
  
  cmd2 = sprintf("cat %s | sort -n | uniq -c > %s", tempFile, fragLenFile)
  system(cmd2)
  
  system(sprintf("rm %s", tempFile))
}

gc()

sessionInfo()
q(save="no")
