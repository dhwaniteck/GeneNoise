library("scran") #for DM function
library("dplyr")
library(org.Sc.sgd.db)

gasch_rna <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/Gasch_raw_data - salt stress/GSE102475_GASCH_NaCl-scRNAseq_NormData.txt", header = TRUE, sep = "\t", dec = ".")

rna_noise_gasch <- as.data.frame(gasch_rna$gene)
colnames(rna_noise_gasch)[1] <- "ORF"
rna_noise_gasch$ORF[!(rna_noise_gasch$ORF %in% keys(org.Sc.sgd.db, keytype = "ORF"))] # all the ERCC are not yeast genes and will be removed from the analysis subsequently

rna_noise_gasch$AVGEXP_stressed <- rowMeans(dplyr::select(gasch_rna, contains("_Stressed_")), na.rm = TRUE)
rna_noise_gasch$SD_stressed <- apply(dplyr::select(gasch_rna, contains("_Stressed_")), 1, sd, na.rm=TRUE)
rna_noise_gasch$CV_stressed <- rna_noise_gasch$SD_stressed/rna_noise_gasch$AVGEXP_stressed
rna_noise_gasch$CV_stressed[which(rowSums(dplyr::select(gasch_rna, contains("_Stressed_")) != 0) < 4)] <- NA   #Calculate CV only for those genes which have non-zero values at least in 4 cells

rna_noise_gasch$AVGEXP_unstressed <- rowMeans(dplyr::select(gasch_rna, contains("_Unstressed_")), na.rm = TRUE)
rna_noise_gasch$SD_unstressed <- apply(dplyr::select(gasch_rna, contains("_Unstressed_")), 1, sd, na.rm=TRUE)
rna_noise_gasch$CV_unstressed <- rna_noise_gasch$SD_unstressed/rna_noise_gasch$AVGEXP_unstressed
rna_noise_gasch$CV_unstressed[which(rowSums(dplyr::select(gasch_rna, contains("_Unstressed_")) != 0) < 4)] <- NA   #Calculate CV only for those genes which have non-zero values at least in 4 cells

#rna_noise_gasch <- rna_noise_gasch[121:6789,]
rna_noise_gasch <- rna_noise_gasch[(rna_noise_gasch$ORF %in% keys(org.Sc.sgd.db, keytype = "ORF")),]

rna_noise_gasch$DM_stressed <- DM(rna_noise_gasch$AVGEXP_stressed, rna_noise_gasch$CV_stressed ^ 2, win.size=100)
rna_noise_gasch$FF_stressed <- (rna_noise_gasch$SD_stressed ^ 2)/rna_noise_gasch$AVGEXP_stressed #Fano factor
#rna_noise_gasch$FF_stressed[which(rowSums(dplyr::select(gasch_rna, contains("_Stressed_")) != 0) < 4)] <- NA
plot(log(rna_noise_gasch$AVGEXP_stressed),log(rna_noise_gasch$CV_stressed^2))
plot(log(rna_noise_gasch$AVGEXP_stressed),rna_noise_gasch$DM_stressed)

rna_noise_gasch$DM_unstressed <- DM(rna_noise_gasch$AVGEXP_unstressed, rna_noise_gasch$CV_unstressed ^ 2, win.size=100)
rna_noise_gasch$FF_unstressed <- (rna_noise_gasch$SD_unstressed ^ 2)/rna_noise_gasch$AVGEXP_unstressed #Fano factor
#rna_noise_gasch$FF_unstressed[which(rowSums(dplyr::select(gasch_rna, contains("A")) != 0) < 4)] <- NA
plot(log(rna_noise_gasch$AVGEXP_unstressed),log(rna_noise_gasch$CV_unstressed^2))
plot(log(rna_noise_gasch$AVGEXP_unstressed),rna_noise_gasch$DM_unstressed)

write.csv(rna_noise_gasch,"/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/RNA_noise_gasch.csv")

# Understanding how the DM function works and what it means to calculate the running median
# test <-  data_frame(cv = sample(1:20,10),
#                        mean = sample(1:20,10))
# plot(test$mean,test$cv)
# test$logcv <- log10(test$cv)
# test$logmean <- log10(test$mean)
# test$dm <- DM(test$mean,test$cv, win.size = 3)
# test$calccv <- test$logcv - test$dm
# o <- order(test$logmean)
# test$medtrend <- runmed(test$logcv[o], k = 3)
# test$medtrend[o] <- test$medtrend
