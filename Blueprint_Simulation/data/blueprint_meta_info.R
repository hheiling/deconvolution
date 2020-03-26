
# The file "blueprint_meta_info.rds" has the meta infomratoin

meta = readRDS("blueprint_meta_info.rds")
dim(meta)
meta[1:2,]

sapply(c("LIBRARY_LAYOUT","LIBRARY_STRATEGY","BIOMATERIAL_TYPE","CELL_TYPE",
         "DONOR_SEX","DONOR_ETHNICITY","TISSUE_TYPE","SPECIMEN_PROCESSING",
         "DATASET_ID","DATASET_TITLE"),
       function(xx) table(meta[,xx], useNA="ifany"))

table(meta$LIBRARY_LAYOUT, meta$DATASET_ID)

sessionInfo()

q(save="no")
