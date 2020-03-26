# Examination of fragment length distributions

library(plotrix)

fragFiles1 = list.files(path = "C:/Users/hheiling/Documents/Longleaf/Blueprint/Fragment_Lengths/EGAD00001002671/",
                        full.names = T)
length(fragFiles1)

# par(mfrow = c(3,3), mar = c(3,4,3,2))
# for(i in 1:9){
#   fDist = read.table(fragFiles1[i], as.is = T)
#   weighted.hist(x = fDist[which(fDist[,1]>1000),2], w = fDist[which(fDist[,1]>1000),1])
# }

df_fragLen = read.table(fragFiles1[1], as.is = T)
colnames(df_fragLen) = c("Freq","Len")
df_fragLen$Weight = df_fragLen$Freq / sum(df_fragLen$Freq)
df_fragLen$SampleID = rep(1, times = nrow(df_fragLen))

for(i in 2:50){
  fDist = read.table(fragFiles1[i], as.is = T)
  fDist[,3] = fDist[,1] / sum(fDist[,1])
  fDist[,4] = rep(i, times = nrow(fDist))
  colnames(fDist) = c("Freq","Len","Weight","SampleID")
  
  df_fragLen = rbind(df_fragLen, fDist)
}

library(ggplot2)
ggplot(data = df_fragLen, mapping = aes(x = Len, weight = Weight, group = SampleID)) + geom_density() +
  coord_cartesian(xlim = c(100, 3*10^4)) + # coord_cartesian(xlim = c(100, 3*10^4))
  ggtitle("FragLen Density Plots for 50 Samples") + xlab("Length of Fragment")

###########################################################################################################

# Combine multiple fragment length files together to create one overall fragment length distribution
# file to be used in the Blueprint model fit

min_len = min(df_fragLen$Len)
min_len
max_len = max(df_fragLen$Len)
max_len

###########################################################################################################
