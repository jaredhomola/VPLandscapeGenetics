###################################################################
######      Implementation of multiple matrix regression     ######
######     with randomization via lgrMMRR of PopGenReport    ######
###################################################################

##### Load required packages #####
library(geosphere)
library(PopGenReport)
library(devtools)
library(VPLandscapeGenetics)
data(VPLandscapeGenetics)


##### Calculate genetic distances #####
##SS
nei.SS <- as.matrix(pairwise_Gst_Nei(genind.SS, linearized = TRUE))
nei.SS <- as.matrix(nei.SS)
dimnames(nei.SS) <- list(1:90, 1:90)

##WF
nei.WF <- as.matrix(pairwise_Gst_Nei(genind.WF, linearized = TRUE))
nei.WF <- as.matrix(nei.WF)
dimnames(nei.WF) <- list(1:87, 1:87)


##### Calculate geographic distances #####
##SS
geoDist.SS <- as.matrix(distm(latLong.SS, fun=distGeo))
geoDist.SS <- geoDist.SS/1000
dim(geoDist.SS) <- c(90, 90)
dimnames(geoDist.SS) <- list(1:90, 1:90)

##WF
geoDist.WF <- as.matrix(distm(latLong.WF, fun=distGeo))
geoDist.WF <- geoDist.WF/1000
dim(geoDist.WF) <- c(87, 87)
dimnames(geoDist.WF) <- list(1:87, 1:87)


##### Calculate distance matricies for environmental metrics #####
##SS
ltRds.SS <- as.numeric(as.matrix(dist(envData.SS$SS_lightroads, method="euclidean")))
secRds.SS <- as.numeric(as.matrix(dist(envData.SS$SS_secondaryroads, method="euclidean")))
prmRds.SS <- as.numeric(as.matrix(dist(envData.SS$SS_primaryroads, method="euclidean")))

dim(ltRds.SS) <- c(90, 90)
dimnames(ltRds.SS) <- list(1:90, 1:90)
dim(secRds.SS) <- c(90, 90)
dimnames(secRds.SS) <- list(1:90, 1:90)
dim(prmRds.SS) <- c(90, 90)
dimnames(prmRds.SS) <- list(1:90, 1:90)

envDist_list.SS <- list(ltRds.SS, secRds.SS, prmRds.SS)

##WF
ltRds.WF <- as.numeric(as.matrix(dist(envData.WF$WF_lightroads, method="euclidean")))
secRds.WF <- as.numeric(as.matrix(dist(envData.WF$WF_secondaryroads, method="euclidean")))
prmRds.WF <- as.numeric(as.matrix(dist(envData.WF$WF_primaryroads, method="euclidean")))

dim(ltRds.WF) <- c(87, 87)
dimnames(ltRds.WF) <- list(1:87, 1:87)
dim(secRds.WF) <- c(87, 87)
dimnames(secRds.WF) <- list(1:87, 1:87)
dim(prmRds.WF) <- c(87, 87)
dimnames(prmRds.WF) <- list(1:87, 1:87)

envDist_list.WF <- list(ltRds.WF, secRds.WF, prmRds.WF)

##### Perform MMRR analyses #####
##SS
mmrr.ss <- lgrMMRR(nei.SS, envDist_list.SS, geoDist.SS, nperm=10000)

##WF
mmrr.WF <- lgrMMRR(nei.WF, envDist_list.WF, geoDist.WF, nperm=10000)
