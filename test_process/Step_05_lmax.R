# Find lmax from fragment lengths output (Step_03_FragLen.R)

file = "set1_50_set2_50" # Should be one of the comboLabels provided in FragLengths_Function_mixture.R

fragSizeFile = sprintf("%s_fraglens.txt",file)

fragDist = read.table(fragSizeFile)

class(fragDist) # data.frame with 2 columns (unlabeled): Freq, Len
head(fragDist)

summary(fragDist[,2])

lmax_option = 600
sum(fragDist[,2] <= lmax_option)
# head(fragDist)
# summary(fragDist)