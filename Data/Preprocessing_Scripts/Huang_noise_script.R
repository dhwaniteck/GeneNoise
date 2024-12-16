library("scran") #for DM function
library("dplyr")
library(org.Sc.sgd.db)
library("AnnotationDbi")

#A: 2hr , B: 16hr , C: 36hr
huang_rna <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/Huang_raw_data - yeast aging/GSE210032_Log2NormCounts_1W-A37B43C45_H2_NoERCC_CSV.csv", header = TRUE, sep = ",", dec = ".")

rna_noise_huang <- as.data.frame(huang_rna$Geneid)
colnames(rna_noise_huang)[1] <- "ORF"

table(rna_noise_huang$ORF %in% keys(org.Sc.sgd.db, keytype = "ORF")) # checking if all rows correspond to genes in the yeast database
rna_noise_huang$ORF[!(rna_noise_huang$ORF %in% keys(org.Sc.sgd.db, keytype = "ORF"))] # This says that only one gene is not in the database, but on investigation that gene is actually a yeast gene, just under a different ORF. So, the gene's (AAD6) ORF (YFL057C) has been updated to YFL056C
AnnotationDbi::select(org.Sc.sgd.db, keys = "YFL057C", columns = "ORF", keytype = "COMMON") #checking the common name of this gene
AnnotationDbi::select(org.Sc.sgd.db, keys = "YFL057C", columns = "GENENAME", keytype = "COMMON") #finding th ORF of this gene
rna_noise_huang$ORF[!(rna_noise_huang$ORF %in% keys(org.Sc.sgd.db, keytype = "ORF"))] <- "YFL056C"
table(rna_noise_huang$ORF %in% keys(org.Sc.sgd.db, keytype = "ORF")) #Now all rows pass the test of being present in the database as yeast genes.

rna_noise_huang$AVGEXP_2hr <- rowMeans(dplyr::select(huang_rna, contains("A")), na.rm = TRUE)
rna_noise_huang$SD_2hr <- apply(dplyr::select(huang_rna, contains("A")), 1, sd, na.rm=TRUE)
rna_noise_huang$CV_2hr <- rna_noise_huang$SD_2hr/rna_noise_huang$AVGEXP_2hr
rna_noise_huang$CV_2hr[which(rowSums(dplyr::select(huang_rna, contains("A")) != 0) < 4)] <- NA   #Calculate CV only for those genes which have non-zero values at least in 4 cells
rna_noise_huang$DM_2hr <- DM(rna_noise_huang$AVGEXP_2hr, rna_noise_huang$CV_2hr ^ 2, win.size=100)
rna_noise_huang$FF_2hr <- (rna_noise_huang$SD_2hr ^ 2)/rna_noise_huang$AVGEXP_2hr #Fano factor
rna_noise_huang$FF_2hr[which(rowSums(dplyr::select(huang_rna, contains("A")) != 0) < 4)] <- NA
plot(log(rna_noise_huang$AVGEXP_2hr),log(rna_noise_huang$CV_2hr^2))
plot(log(rna_noise_huang$AVGEXP_2hr[rna_noise_huang$AVGEXP_2hr>0.1]),rna_noise_huang$DM_2hr[rna_noise_huang$AVGEXP_2hr>0.1])

rna_noise_huang$AVGEXP_16hr <- rowMeans(dplyr::select(huang_rna, contains("B")), na.rm = TRUE)
rna_noise_huang$SD_16hr <- apply(dplyr::select(huang_rna, contains("B")), 1, sd, na.rm=TRUE)
rna_noise_huang$CV_16hr <- rna_noise_huang$SD_16hr/rna_noise_huang$AVGEXP_16hr
rna_noise_huang$CV_16hr[which(rowSums(dplyr::select(huang_rna, contains("B")) != 0) < 4)] <- NA
rna_noise_huang$DM_16hr <- DM(rna_noise_huang$AVGEXP_16hr, rna_noise_huang$CV_16hr ^ 2, win.size=100)
rna_noise_huang$FF_16hr <- (rna_noise_huang$SD_16hr ^ 2)/rna_noise_huang$AVGEXP_16hr #Fano factor
rna_noise_huang$FF_16hr[which(rowSums(dplyr::select(huang_rna, contains("B")) != 0) < 4)] <- NA
plot(log(rna_noise_huang$AVGEXP_16hr),log(rna_noise_huang$CV_16hr^2))
plot(log(rna_noise_huang$AVGEXP_16hr[rna_noise_huang$AVGEXP_16hr>0.1]),rna_noise_huang$DM_16hr[rna_noise_huang$AVGEXP_16hr>0.1])

rna_noise_huang$AVGEXP_36hr <- rowMeans(dplyr::select(huang_rna, contains("C")), na.rm = TRUE)
rna_noise_huang$SD_36hr <- apply(dplyr::select(huang_rna, contains("C")), 1, sd, na.rm=TRUE)
rna_noise_huang$CV_36hr <- rna_noise_huang$SD_36hr/rna_noise_huang$AVGEXP_36hr
rna_noise_huang$CV_36hr[which(rowSums(dplyr::select(huang_rna, contains("C")) != 0) < 4)] <- NA
rna_noise_huang$DM_36hr <- DM(rna_noise_huang$AVGEXP_36hr, rna_noise_huang$CV_36hr ^ 2, win.size=100)
rna_noise_huang$FF_36hr <- (rna_noise_huang$SD_36hr ^ 2)/rna_noise_huang$AVGEXP_36hr #Fano factor
rna_noise_huang$FF_36hr[which(rowSums(dplyr::select(huang_rna, contains("C")) != 0) < 4)] <- NA
plot(log(rna_noise_huang$AVGEXP_36hr[rna_noise_huang$AVGEXP_36hr>0.1]),log(rna_noise_huang$CV_36hr[rna_noise_huang$AVGEXP_36hr>0.1]^2))
plot(log(rna_noise_huang$AVGEXP_36hr[rna_noise_huang$AVGEXP_36hr>0.1]),rna_noise_huang$DM_36hr[rna_noise_huang$AVGEXP_36hr>0.1])
abline(1.5,-1)

plot(rna_noise_huang$DM_36hr,rna_noise_huang$DM_16hr)
plot(rna_noise_huang$DM_36hr,rna_noise_huang$DM_2hr)
cor.test(rna_noise_huang$DM_36hr,rna_noise_huang$DM_2hr)

boxplot(rna_noise_huang[,c('FF_2hr','FF_16hr','FF_36hr')])
hist(rna_noise_huang$FF_2hr,breaks = 100)
plot(rna_noise_huang$AVGEXP_2hr,rna_noise_huang$FF_2hr)
table(rna_noise_huang$FF_2hr>1.5)

write.csv(rna_noise_huang,"/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/RNA_noise_huang.csv")
