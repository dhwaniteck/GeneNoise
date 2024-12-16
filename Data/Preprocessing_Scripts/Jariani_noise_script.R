library("scran") #for DM function
library("dplyr")
library(org.Sc.sgd.db)

# No. of cells per condition
# glucose 6hr = 1696
# glucose 12hr = 1097
# lag 1hr = 1173
# lag 3hr = 1470
# glucose + maltose = 687
jariani_rna_glu6hr <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/Jariani_raw_data - GSE144820_RAW - glucose-maltose lag/GSM4297054_processed_counts_glu_6h.csv", header = TRUE, sep = ",", dec = ".")
jariani_rna_glu12hr <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/Jariani_raw_data - GSE144820_RAW - glucose-maltose lag/GSM4297055_processed_counts_glu_12h.csv", header = TRUE, sep = ",", dec = ".")
jariani_rna_lag1hr <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/Jariani_raw_data - GSE144820_RAW - glucose-maltose lag/GSM4297056_processed_counts_inlag_1h.csv", header = TRUE, sep = ",", dec = ".")
jariani_rna_lag3hr <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/Jariani_raw_data - GSE144820_RAW - glucose-maltose lag/GSM4297057_processed_counts_inlag_3h.csv", header = TRUE, sep = ",", dec = ".")
jariani_rna_glumal <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/Jariani_raw_data - GSE144820_RAW - glucose-maltose lag/GSM4297058_processed_counts_mix_glu_mal.csv", header = TRUE, sep = ",", dec = ".")

rna_noise_jariani <- as.data.frame(jariani_rna_glu6hr$X)
colnames(rna_noise_jariani)[1] <- "Gene"
rna_noise_jariani$ORF <- mapIds(org.Sc.sgd.db, 
                                keys = rna_noise_jariani$Gene, 
                                column = "ORF", 
                                keytype = "GENENAME", 
                                multiVals = "first")
rna_noise_jariani$ORF[is.na(rna_noise_jariani$ORF)] <- rna_noise_jariani$Gene[is.na(rna_noise_jariani$ORF)]
rna_noise_jariani$Gene <- NULL

rna_noise_jariani$AVGEXP_glu6hr <- rowMeans(jariani_rna_glu6hr[,-1], na.rm = TRUE)
rna_noise_jariani$SD_glu6hr <- apply(jariani_rna_glu6hr[,-1], 1, sd, na.rm=TRUE)
rna_noise_jariani$CV_glu6hr <- rna_noise_jariani$SD_glu6hr/rna_noise_jariani$AVGEXP_glu6hr
rna_noise_jariani$CV_glu6hr[which(rowSums(jariani_rna_glu6hr[,-1]) != 0) < 100] <- NA   #Calculate CV only for those genes which have non-zero values at least in 4 cells

rna_noise_jariani$AVGEXP_glu12hr <- rowMeans(jariani_rna_glu12hr[,-1], na.rm = TRUE)
rna_noise_jariani$SD_glu12hr <- apply(jariani_rna_glu12hr[,-1], 1, sd, na.rm=TRUE)
rna_noise_jariani$CV_glu12hr <- rna_noise_jariani$SD_glu12hr/rna_noise_jariani$AVGEXP_glu12hr
rna_noise_jariani$CV_glu12hr[which(rowSums(jariani_rna_glu12hr[,-1]) != 0) < 100] <- NA   #Calculate CV only for those genes which have non-zero values at least in 4 cells

rna_noise_jariani$AVGEXP_lag1hr <- rowMeans(jariani_rna_lag1hr[,-1], na.rm = TRUE)
rna_noise_jariani$SD_lag1hr <- apply(jariani_rna_lag1hr[,-1], 1, sd, na.rm=TRUE)
rna_noise_jariani$CV_lag1hr <- rna_noise_jariani$SD_lag1hr/rna_noise_jariani$AVGEXP_lag1hr
rna_noise_jariani$CV_lag1hr[which(rowSums(jariani_rna_lag1hr[,-1]) != 0) < 100] <- NA   #Calculate CV only for those genes which have non-zero values at least in 4 cells

rna_noise_jariani$AVGEXP_lag3hr <- rowMeans(jariani_rna_lag3hr[,-1], na.rm = TRUE)
rna_noise_jariani$SD_lag3hr <- apply(jariani_rna_lag3hr[,-1], 1, sd, na.rm=TRUE)
rna_noise_jariani$CV_lag3hr <- rna_noise_jariani$SD_lag3hr/rna_noise_jariani$AVGEXP_lag3hr
rna_noise_jariani$CV_lag3hr[which(rowSums(jariani_rna_lag3hr[,-1]) != 0) < 100] <- NA   #Calculate CV only for those genes which have non-zero values at least in 4 cells

rna_noise_jariani$AVGEXP_glumal <- rowMeans(jariani_rna_glumal[,-1], na.rm = TRUE)
rna_noise_jariani$SD_glumal <- apply(jariani_rna_glumal[,-1], 1, sd, na.rm=TRUE)
rna_noise_jariani$CV_glumal <- rna_noise_jariani$SD_glumal/rna_noise_jariani$AVGEXP_glumal
rna_noise_jariani$CV_glumal[which(rowSums(jariani_rna_glumal[,-1]) != 0) < 100] <- NA   #Calculate CV only for those genes which have non-zero values at least in 4 cells

rna_noise_jariani <- rna_noise_jariani[(rna_noise_jariani$ORF %in% keys(org.Sc.sgd.db, keytype = "ORF")),] #remove genes which are not found in the yeast database

rna_noise_jariani$DM_glu6hr <- DM(rna_noise_jariani$AVGEXP_glu6hr, rna_noise_jariani$CV_glu6hr ^ 2, win.size=100)
rna_noise_jariani$FF_glu6hr <- (rna_noise_jariani$SD_glu6hr ^ 2)/rna_noise_jariani$AVGEXP_glu6hr #Fano factor
rna_noise_jariani$FF_glu6hr[which(rowSums(jariani_rna_glu6hr[,-1]) != 0) < 100] <- NA
plot(log(rna_noise_jariani$AVGEXP_glu6hr),log(rna_noise_jariani$CV_glu6hr^2))
plot(log(rna_noise_jariani$AVGEXP_glu6hr),rna_noise_jariani$DM_glu6hr)

rna_noise_jariani$DM_glu12hr <- DM(rna_noise_jariani$AVGEXP_glu12hr, rna_noise_jariani$CV_glu12hr ^ 2, win.size=100)
rna_noise_jariani$FF_glu12hr <- (rna_noise_jariani$SD_glu12hr ^ 2)/rna_noise_jariani$AVGEXP_glu12hr #Fano factor
rna_noise_jariani$FF_glu12hr[which(rowSums(jariani_rna_glu12hr[,-1]) != 0) < 100] <- NA
plot(log(rna_noise_jariani$AVGEXP_glu12hr),log(rna_noise_jariani$CV_glu12hr^2))
plot(log(rna_noise_jariani$AVGEXP_glu12hr),rna_noise_jariani$DM_glu12hr)

rna_noise_jariani$DM_lag1hr <- DM(rna_noise_jariani$AVGEXP_lag1hr, rna_noise_jariani$CV_lag1hr ^ 2, win.size=100)
rna_noise_jariani$FF_lag1hr <- (rna_noise_jariani$SD_lag1hr ^ 2)/rna_noise_jariani$AVGEXP_lag1hr #Fano factor
rna_noise_jariani$FF_lag1hr[which(rowSums(jariani_rna_lag1hr[,-1]) != 0) < 100] <- NA
plot(log(rna_noise_jariani$AVGEXP_lag1hr),log(rna_noise_jariani$CV_lag1hr^2))
plot(log(rna_noise_jariani$AVGEXP_lag1hr),rna_noise_jariani$DM_lag1hr)

rna_noise_jariani$DM_lag3hr <- DM(rna_noise_jariani$AVGEXP_lag3hr, rna_noise_jariani$CV_lag3hr ^ 2, win.size=100)
rna_noise_jariani$FF_lag3hr <- (rna_noise_jariani$SD_lag3hr ^ 2)/rna_noise_jariani$AVGEXP_lag3hr #Fano factor
rna_noise_jariani$FF_lag3hr[which(rowSums(jariani_rna_lag3hr[,-1]) != 0) < 100] <- NA
plot(log(rna_noise_jariani$AVGEXP_lag3hr),log(rna_noise_jariani$CV_lag3hr^2))
plot(log(rna_noise_jariani$AVGEXP_lag3hr),rna_noise_jariani$DM_lag3hr)

rna_noise_jariani$DM_glumal <- DM(rna_noise_jariani$AVGEXP_glumal, rna_noise_jariani$CV_glumal ^ 2, win.size=100)
rna_noise_jariani$FF_glumal <- (rna_noise_jariani$SD_glumal ^ 2)/rna_noise_jariani$AVGEXP_glumal #Fano factor
rna_noise_jariani$FF_glumal[which(rowSums(jariani_rna_glumal[,-1]) != 0) < 100] <- NA
plot(log(rna_noise_jariani$AVGEXP_glumal),log(rna_noise_jariani$CV_glumal^2))
plot(log(rna_noise_jariani$AVGEXP_glumal),rna_noise_jariani$DM_glumal)


boxplot(rna_noise_jariani[,c('FF_glu6hr','FF_glu12hr','FF_lag1hr','FF_lag3hr','FF_glumal')])
boxplot(rna_noise_jariani[,c('DM_glu6hr','DM_glu12hr','DM_lag1hr','DM_lag3hr','DM_glumal')])
hist(rna_noise_jariani$FF_glu6hr, breaks = 30)
table(rna_noise_jariani$FF_glu6hr < 1.5 & rna_noise_jariani$FF_glu6hr > 0.5)
table(rna_noise_jariani$FF_glu12hr < 1.5 & rna_noise_jariani$FF_glu12hr > 0.5)

write.csv(rna_noise_jariani,"/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/RNA_noise_jariani.csv")

library(ggfortify)
res.pca <- prcomp(jariani_rna_glu6hr[,-1], scale = T)
autoplot(res.pca)
