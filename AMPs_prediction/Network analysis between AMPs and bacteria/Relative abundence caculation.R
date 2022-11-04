# integration part
name <- read.csv("name.csv", header = F) # name.csv contain re_ab.csv for each line
n <- nrow(name)
na <- t(name)
result <- c(0)
for(i in 1:n){
data <- read.csv(as.character(name[i,]), header = F, sep = '\t')
data <- data[-nrow(data),]
result <- cbind(result, data[,3])
}
result <- as.data.frame(result)
result$result <- data[,1]
na <- cbind('Sample', na)
colnames(result) <- na
# 0.05 part
nr <- nrow(result)
rn <- c(0)
limit <- round(0.05*(ncol(result) - 1)) # at lest present in 5% samples
for(j in 1:nr){
l <- length(which((result[j,] != 0)))
if(l <= limit){
  rn <- c(rn, j)
}
}
rn <- rn[-1]
result <- result[-rn,]
col_n <- ncol(result)
for(ii in 2:col_n){
result[,ii] <- result[,ii]/sum(result[,ii])
}
write.csv(result, "amp0.05.csv", row.names  = F)

# caculates correlation
amp <- read.csv("amp0.05.csv", header = T, sep = ',', row.names = 1)
spe <- read.csv("merged_abundance_table_genus_cut0.05.tsv", header = T, sep = '\t', row.names = 1)
spe <- spe/100 # in file spe, each col sum = 100.
library("WGCNA")
library("multtest")
amp <- t(amp)
spe <- t(spe)
corr <- corAndPvalue(amp, spe, alternative = c("two.sided"), method = "spearman")
# p value adjust part
mtadj <- mt.rawp2adjp(unlist(corr$p), proc = "BH")
adpcor <- mtadj$adjp[order(mtadj$index), 2]
occor.p <- matrix(adpcor, dim(amp)[2])
pvalue <- as.data.frame(occor.p)
# remove the positive correlation part
cor <- as.data.frame(corr$cor)
rownames(pvalue) <- row.names(cor)
colnames(pvalue) <- colnames(cor)
cor[cor > 0] <- 0
cor[pvalue > 0.05] <- 0
write.csv(pvalue, "Pvalues-spearman-genus_metasub.csv",row.names = T)
write.csv(cor, "Correlation-spearman-genus_metasub.csv",row.names = T)
