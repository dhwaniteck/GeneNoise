library("scran") #for DM function

acar_rna <- read.delim("/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/Acar_raw_data/Acar_rna_noise.csv", header = TRUE, sep = ",", dec = ".")

acar_rna$DMSO_DM <- DM(acar_rna$AVG_EXP_DMSO, acar_rna$CV_DMSO ^ 2, win.size=100)
plot(log(acar_rna$AVG_EXP_DMSO),acar_rna$CV_DMSO)
plot(log(acar_rna$AVG_EXP_DMSO),acar_rna$DMSO_DM)

acar_rna$Gua_DM <- DM(acar_rna$AVG_EXP_Gua, acar_rna$CV_Gua ^ 2, win.size=100)
plot(log(acar_rna$AVG_EXP_Gua),acar_rna$CV_Gua)
plot(log(acar_rna$AVG_EXP_Gua),acar_rna$Gua_DM)

acar_rna$MPA_DM <- DM(acar_rna$AVG_EXP_MPA, acar_rna$CV_MPA ^ 2, win.size=100)
plot(log(acar_rna$AVG_EXP_MPA),acar_rna$CV_MPA)
plot(log(acar_rna$AVG_EXP_MPA),acar_rna$MPA_DM)

acar_rna$GuaMPA_DM <- DM(acar_rna$AVG_EXP_GuaMPA, acar_rna$CV_GuaMPA ^ 2, win.size=100)
plot(log(acar_rna$AVG_EXP_GuaMPA),acar_rna$CV_GuaMPA)
plot(log(acar_rna$AVG_EXP_GuaMPA),acar_rna$GuaMPA_DM)

acar_rna$Unnamed..0 <- NULL

write.csv(acar_rna, "/Users/dhwani/Documents/Coursework/SEM3/ML for functional genomics - Prof. David Knowles/Project/Data/Acar_rna_noise.csv", row.names = FALSE)

#checking how noise varies across media
plot(acar_rna$DMSO_DM,acar_rna$Gua_DM)
plot(acar_rna$DMSO_DM,acar_rna$MPA_DM)
cor.test(acar_rna$DMSO_DM,acar_rna$MPA_DM)
plot(acar_rna$DMSO_DM,acar_rna$GuaMPA_DM)
cor.test(acar_rna$DMSO_DM,acar_rna$GuaMPA_DM)


cor.test(acar_rna$DMSO_DM,acar_rna$Gua_DM)
lm(acar_rna$DMSO_DM~acar_rna$Gua_DM)
