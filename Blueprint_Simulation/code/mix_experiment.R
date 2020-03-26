
# Experimenting with properties of pure sample files in order to brainstorm how best to 
# combine them to make mixture samples

files1 = list.files(path = "C:/Users/hheiling/Documents/Longleaf/Blueprint/Fit_Samples/EGAD00001002671/",
                   pattern = "trec_exon_set.txt", full.names = T)

length(files1)

files2 = list.files(path = "C:/Users/hheiling/Documents/Longleaf/Blueprint/Fit_Samples/EGAD00001002674/",
                   pattern = "trec_exon_set.txt", full.names = T)

length(files2)

# One sample from two of the cell types
f1 = read.table(files1[1])
f2 = read.table(files2[1])

head(f1)
head(f2)

tail(f1)
tail(f2)

# Proof that both all.equal and identical are only true if both the values and 
# the order of the values are the same between two options
all.equal(c("a","b","c"),c("b","c","a"))
identical(c("a","b","c"),c("b","c","a"))

# Proof that ordering of the rows of the pure count files are the same for all pure count files
all.equal(rownames(f1),rownames(f2))
identical(rownames(f1),rownames(f2))

# Distribution of total counts (paired data, files1)

tot_cts = numeric(length(files1))
for(i in 1:length(files1)){
  ff = read.table(files1[i], as.is = T)
  tot_cts[i] = sum(ff)
}

# total counts approx. normally distributed
hist(tot_cts)
mean(tot_cts)

# Average total counts around 12.2 million
mean(tot_cts) / 10^6
# SD of total counts around 2.7 million
sd(tot_cts) / 10^6
