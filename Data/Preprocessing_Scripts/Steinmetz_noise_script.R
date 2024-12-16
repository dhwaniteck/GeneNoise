library("scran") #for DM function
library("dplyr")
library(org.Sc.sgd.db)

# no. of cells for each condition
# BY4741 in YPD (labelled as YPD): 127
# YJM789 in YPD (labelled as YJM): 48 # YJM is a less- commonly used haploid S. cerevisiae strain that was isolated from an HIV patient

#Load raw scRNA-seq reads
steinmetz_rna_YPD  <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/Steinmetz_raw_data/GSE122392_normExpression_YPD_NatMicro.tab", header = TRUE, sep = "\t", dec = ".")
steinmetz_rna_YJM  <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/Steinmetz_raw_data/GSE122392_normExpression_YJM_NatMicro.tab", header = TRUE, sep = "\t", dec = ".")

rna_noise_steinmetz <- as.data.frame(steinmetz_rna_YPD$geneName)
colnames(rna_noise_steinmetz)[1] <- "ORF"
table(rna_noise_steinmetz$ORF %in% keys(org.Sc.sgd.db, keytype = "ORF")) #not all the rows in the dataset are genes, some are non-coding RNAs - we will remove them from the dataset post processing. There are also rows with two ORFs written separate by a ",". I am not sure why that is but these have also been removed due to ambiguity.
rna_noise_steinmetz$ORF[!(rna_noise_steinmetz$ORF %in% keys(org.Sc.sgd.db, keytype = "ORF"))]

#Calculate mean expression level, std dev in expression and coefficient of variation for all rows (genes + ncRNAs)
rna_noise_steinmetz$AVGEXP_YPD <- rowMeans(steinmetz_rna_YPD[,-c(1,2)], na.rm = TRUE)
rna_noise_steinmetz$SD_YPD <- apply(steinmetz_rna_YPD[,-c(1,2)], 1, sd, na.rm=TRUE)
rna_noise_steinmetz$CV_YPD <- rna_noise_steinmetz$SD_YPD/rna_noise_steinmetz$AVGEXP_YPD
rna_noise_steinmetz$CV_YPD[which(rowSums(steinmetz_rna_YPD[,-c(1,2)]) != 0) < 4] <- NA   #Calculate CV only for those genes which have non-zero values at least in 4 cells

rna_noise_steinmetz$AVGEXP_YJM <- rowMeans(steinmetz_rna_YJM[,-c(1,2)], na.rm = TRUE)
rna_noise_steinmetz$SD_YJM <- apply(steinmetz_rna_YJM[,-c(1,2)], 1, sd, na.rm=TRUE)
rna_noise_steinmetz$CV_YJM <- rna_noise_steinmetz$SD_YJM/rna_noise_steinmetz$AVGEXP_YJM
rna_noise_steinmetz$CV_YJM[which(rowSums(steinmetz_rna_YJM[,-c(1,2)]) != 0) < 4] <- NA   #Calculate CV only for those genes which have non-zero values at least in 4 cells

#filter out ncRNAs
rna_noise_steinmetz <- rna_noise_steinmetz[(rna_noise_steinmetz$ORF %in% keys(org.Sc.sgd.db, keytype = "ORF")),]

#Calculate distance from median (DM ~ Noise) and fano factor (FF)
rna_noise_steinmetz$DM_YPD <- DM(rna_noise_steinmetz$AVGEXP_YPD, rna_noise_steinmetz$CV_YPD ^ 2, win.size=100)
rna_noise_steinmetz$FF_YPD <- (rna_noise_steinmetz$SD_YPD ^ 2)/rna_noise_steinmetz$AVGEXP_YPD #Fano factor

rna_noise_steinmetz$DM_YJM <- DM(rna_noise_steinmetz$AVGEXP_YJM, rna_noise_steinmetz$CV_YJM ^ 2, win.size=100)
rna_noise_steinmetz$FF_YJM <- (rna_noise_steinmetz$SD_YJM ^ 2)/rna_noise_steinmetz$AVGEXP_YJM #Fano factor

plot(log(rna_noise_steinmetz$AVGEXP_YPD),log(rna_noise_steinmetz$CV_YPD^2))
plot(log(rna_noise_steinmetz$AVGEXP_YPD),rna_noise_steinmetz$DM_YPD)
plot(log(rna_noise_steinmetz$AVGEXP_YJM),log(rna_noise_steinmetz$CV_YJM^2))
plot(log(rna_noise_steinmetz$AVGEXP_YJM),rna_noise_steinmetz$DM_YJM)

boxplot(rna_noise_steinmetz[,c('FF_YPD','FF_YJM')])
boxplot(rna_noise_steinmetz[,c('DM_YPD','DM_YJM')])

plot(log(rna_noise_steinmetz$AVGEXP_YPD),rna_noise_steinmetz$FF_YPD)

plot(rna_noise_steinmetz$DM_YPD,rna_noise_steinmetz$DM_YJM)
abline(a=0,b=1)
cor.test(rna_noise_steinmetz$DM_YPD,rna_noise_steinmetz$DM_YJM)

write.csv(rna_noise_steinmetz,"/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/RNA_noise_steinmetz.csv")

library(ggfortify)
res.pca <- prcomp(steinmetz_rna_YPD[,c(-1,-2)], scale = T)
autoplot(res.pca)

res.pca

library(umap)
rna_umap <- umap(steinmetz_rna_YPD[,c(-1,-2)])
plot(x = rna_umap$layout[,1],
     y = rna_umap$layout[,2])
