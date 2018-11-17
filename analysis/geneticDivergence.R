########################################################
#####     Estimation of genetic divergence using   #####
#####    Gst (Nei 1973) and G''st (Hedrick 2005)   #####
########################################################

##### Load required packages #####
library(adegenet)
library(mmod)
library(VPLandscapeGenetics)
data("VPLandscapeGenetics")

##### Esimate population pairwise Gst (Nei) and G''st (Hedrick) #####
hedrick.SS <- pairwise_Gst_Hedrick(genind.SS)
nei.SS <- pairwise_Gst_Nei(genind.SS)

hedrick.WF <- pairwise_Gst_Hedrick(genind.WF)
nei.WF <- pairwise_Gst_Nei(genind.WF)

##### Esimate global Gst (Nei) and G''st (Hedrick) #####
diffStats.SS <- diff_stats(genind.SS, phi_st=FALSE)

diffStats.WF <- diff_stats(genind.WF, phi_st=FALSE)


##### Assemble full matricies for publication ######
## Spotted salamanders
fullMatrix.SS <- matrix(NA, nrow = 90, ncol = 90)
fullMatrix.SS[upper.tri(fullMatrix.SS)] <- as.matrix(hedrick.SS)[upper.tri(hedrick.SS)]
fullMatrix.SS[lower.tri(fullMatrix.SS)] <- as.matrix(nei.SS)[lower.tri(nei.SS)]
rownames(fullMatrix.SS) <- attributes(nei.SS)$Labels
colnames(fullMatrix.SS) <- attributes(nei.SS)$Labels

write.csv(fullMatrix.SS, "differentiationMatrix.SS.csv")

## Wood frogs
fullMatrix.WF <- matrix(NA, nrow = 87, ncol = 87)
fullMatrix.WF[upper.tri(fullMatrix.WF)] <- as.matrix(hedrick.WF)[upper.tri(hedrick.WF)]
fullMatrix.WF[lower.tri(fullMatrix.WF)] <- as.matrix(nei.WF)[lower.tri(nei.WF)]
rownames(fullMatrix.WF) <- attributes(nei.WF)$Labels
colnames(fullMatrix.WF) <- attributes(nei.WF)$Labels
write.csv(fullMatrix.WF, "differentiationMatrixWF.csv")

