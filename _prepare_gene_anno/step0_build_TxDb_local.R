
library("GenomicFeatures")

gtfFile = "gencode.v15.annotation.gtf.gz"

path = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_15/"

txdb = makeTxDbFromGFF(file=gtfFile, format="gtf",
  dataSource=paste(path, gtfFile, sep=""),
  organism="Homo sapiens")

seqlevels(txdb)
columns(txdb)
keytypes(txdb)

genes = exonsBy(txdb, by="gene")
saveRDS(genes, file = "exon_by_genes_gencode.v15.GRCh37.rds")

sessionInfo()

q(save="no")
