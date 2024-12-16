library("scran") #for DM function
library("dplyr")
library(org.Sc.sgd.db)
library(Matrix) # to read .dge file
library("AnnotationDbi")

#Acar lab data - 4 samples (DMSO n=233, Guanine n=258, MPA n=85, Guanine MPA n=268 cells)

# Read the .mtx file
acar_rna <- readMM("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/Acar_raw_data/GSE165686_yeastdropseq_dge.mtx")
acar_rna <- as.matrix(acar_rna)
genes <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/Acar_raw_data/GSE165686_features.tsv", header = TRUE, sep = "\t", dec = ".")
genes$x[!(genes$x %in% keys(org.Sc.sgd.db, keytype = "ORF"))]
AnnotationDbi::select(org.Sc.sgd.db, keys = "YFL057C", columns = "ORF", keytype = "COMMON") #checking the common name of this gene
AnnotationDbi::select(org.Sc.sgd.db, keys = "YFL057C", columns = "GENENAME", keytype = "COMMON") #finding th ORF of this gene
AnnotationDbi::select(org.Sc.sgd.db, keys = "YIL168W", columns = "ORF", keytype = "COMMON") #checking the common name of this gene
AnnotationDbi::select(org.Sc.sgd.db, keys = "YIL168W", columns = "GENENAME", keytype = "COMMON") #finding th ORF of this gene
genes$x[!(genes$x %in% keys(org.Sc.sgd.db, keytype = "ORF"))] <- c("YFL056C","YIL167W")
genes$x[!(genes$x %in% keys(org.Sc.sgd.db, keytype = "ORF"))] # Now all genes are confirmed to have the correct ORFs

rna_noise_acar <- as.data.frame(genes$x)
colnames(rna_noise_acar)[1] <- "ORF"

rna_noise_acar$AVGEXP_DMSO <- rowMeans(acar_rna[,1:233], na.rm = TRUE)
rna_noise_acar$SD_DMSO <- apply(acar_rna[,1:233], 1, sd, na.rm=TRUE)
rna_noise_acar$CV_DMSO <- rna_noise_acar$SD_DMSO/rna_noise_acar$AVGEXP_DMSO
rna_noise_acar$CV_DMSO[which(rowSums(acar_rna[,1:233] != 0) < 4)] <- NA   #Calculate CV only for those genes which have non-zero values at least in 4 cells
rna_noise_acar$DM_DMSO <- DM(rna_noise_acar$AVGEXP_DMSO, rna_noise_acar$CV_DMSO ^ 2, win.size=100)
rna_noise_acar$FF_DMSO <- (rna_noise_acar$SD_DMSO ^ 2)/rna_noise_acar$AVGEXP_DMSO #Fano factor
rna_noise_acar$FF_DMSO[which(rowSums(acar_rna[,1:233] != 0) < 4)] <- NA
plot(log(rna_noise_acar$AVGEXP_DMSO),log(rna_noise_acar$CV_DMSO^2))
plot(log(rna_noise_acar$AVGEXP_DMSO),rna_noise_acar$DM_DMSO)

rna_noise_acar$AVGEXP_gua <- rowMeans(acar_rna[,234:491], na.rm = TRUE)
rna_noise_acar$SD_gua <- apply(acar_rna[,234:491], 1, sd, na.rm=TRUE)
rna_noise_acar$CV_gua <- rna_noise_acar$SD_gua/rna_noise_acar$AVGEXP_gua
rna_noise_acar$CV_gua[which(rowSums(acar_rna[,234:491] != 0) < 4)] <- NA   #Calculate CV only for those genes which have non-zero values at least in 4 cells
rna_noise_acar$DM_gua <- DM(rna_noise_acar$AVGEXP_gua, rna_noise_acar$CV_gua ^ 2, win.size=100)
rna_noise_acar$FF_gua <- (rna_noise_acar$SD_gua ^ 2)/rna_noise_acar$AVGEXP_gua #Fano factor
rna_noise_acar$FF_gua[which(rowSums(acar_rna[,234:491] != 0) < 4)] <- NA
plot(log(rna_noise_acar$AVGEXP_gua),log(rna_noise_acar$CV_gua^2))
plot(log(rna_noise_acar$AVGEXP_gua),rna_noise_acar$DM_gua)

rna_noise_acar$AVGEXP_MPA <- rowMeans(acar_rna[,492:576], na.rm = TRUE)
rna_noise_acar$SD_MPA <- apply(acar_rna[,492:576], 1, sd, na.rm=TRUE)
rna_noise_acar$CV_MPA <- rna_noise_acar$SD_MPA/rna_noise_acar$AVGEXP_MPA
rna_noise_acar$CV_MPA[which(rowSums(acar_rna[,492:576] != 0) < 4)] <- NA   #Calculate CV only for those genes which have non-zero values at least in 4 cells
rna_noise_acar$DM_MPA <- DM(rna_noise_acar$AVGEXP_MPA, rna_noise_acar$CV_MPA ^ 2, win.size=100)
rna_noise_acar$FF_MPA <- (rna_noise_acar$SD_MPA ^ 2)/rna_noise_acar$AVGEXP_MPA #Fano factor
rna_noise_acar$FF_MPA[which(rowSums(acar_rna[,492:576] != 0) < 4)] <- NA
plot(log(rna_noise_acar$AVGEXP_MPA),log(rna_noise_acar$CV_MPA^2))
plot(log(rna_noise_acar$AVGEXP_MPA),rna_noise_acar$DM_MPA)

rna_noise_acar$AVGEXP_guaMPA <- rowMeans(acar_rna[,577:844], na.rm = TRUE)
rna_noise_acar$SD_guaMPA <- apply(acar_rna[,577:844], 1, sd, na.rm=TRUE)
rna_noise_acar$CV_guaMPA <- rna_noise_acar$SD_guaMPA/rna_noise_acar$AVGEXP_guaMPA
rna_noise_acar$CV_guaMPA[which(rowSums(acar_rna[,577:844] != 0) < 4)] <- NA   #Calculate CV only for those genes which have non-zero values at least in 4 cells
rna_noise_acar$DM_guaMPA <- DM(rna_noise_acar$AVGEXP_guaMPA, rna_noise_acar$CV_guaMPA ^ 2, win.size=100)
rna_noise_acar$FF_guaMPA <- (rna_noise_acar$SD_guaMPA ^ 2)/rna_noise_acar$AVGEXP_guaMPA #Fano factor
rna_noise_acar$FF_guaMPA[which(rowSums(acar_rna[,577:844] != 0) < 4)] <- NA
plot(log(rna_noise_acar$AVGEXP_guaMPA),log(rna_noise_acar$CV_guaMPA^2))
plot(log(rna_noise_acar$AVGEXP_guaMPA),rna_noise_acar$DM_guaMPA)

write.csv(rna_noise_acar, "/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/RNA_noise_acar.csv", row.names = FALSE)


# Initial script when partial processing was done in python --------
# 
# acar_rna <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/Acar_raw_data/Acar_rna_noise.csv", header = TRUE, sep = ",", dec = ".")
# 
# acar_rna$DMSO_DM <- DM(acar_rna$AVG_EXP_DMSO, acar_rna$CV_DMSO ^ 2, win.size=100)
# plot(log(acar_rna$AVG_EXP_DMSO),acar_rna$CV_DMSO)
# plot(log(acar_rna$AVG_EXP_DMSO),acar_rna$DMSO_DM)
# 
# acar_rna$Gua_DM <- DM(acar_rna$AVG_EXP_Gua, acar_rna$CV_Gua ^ 2, win.size=100)
# plot(log(acar_rna$AVG_EXP_Gua),acar_rna$CV_Gua)
# plot(log(acar_rna$AVG_EXP_Gua),acar_rna$Gua_DM)
# 
# acar_rna$MPA_DM <- DM(acar_rna$AVG_EXP_MPA, acar_rna$CV_MPA ^ 2, win.size=100)
# plot(log(acar_rna$AVG_EXP_MPA),acar_rna$CV_MPA)
# plot(log(acar_rna$AVG_EXP_MPA),acar_rna$MPA_DM)
# 
# acar_rna$GuaMPA_DM <- DM(acar_rna$AVG_EXP_GuaMPA, acar_rna$CV_GuaMPA ^ 2, win.size=100)
# plot(log(acar_rna$AVG_EXP_GuaMPA),acar_rna$CV_GuaMPA)
# plot(log(acar_rna$AVG_EXP_GuaMPA),acar_rna$GuaMPA_DM)
# 
# acar_rna$Unnamed..0 <- NULL
# 
# write.csv(acar_rna, "/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/Acar_rna_noise.csv", row.names = FALSE)
# 
# #checking how noise varies across media
# plot(acar_rna$DMSO_DM,acar_rna$Gua_DM)
# plot(acar_rna$DMSO_DM,acar_rna$MPA_DM)
# cor.test(acar_rna$DMSO_DM,acar_rna$MPA_DM)
# plot(acar_rna$DMSO_DM,acar_rna$GuaMPA_DM)
# cor.test(acar_rna$DMSO_DM,acar_rna$GuaMPA_DM)
# 
# 
# cor.test(acar_rna$DMSO_DM,acar_rna$Gua_DM)
# lm(acar_rna$DMSO_DM~acar_rna$Gua_DM)
