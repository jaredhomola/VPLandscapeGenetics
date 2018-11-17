########################################################
####      Perform Mantel tests for both species.    ####
########################################################
########################################################
#####  Genetic distances included Gst (Nei 1973)   #####
#####           and G''st (Hedrick 2005)           #####
########################################################
#####  Geographic distances were evaluated with    #####
#####       and without log transformation         #####
#####         (i.e., Rousset 1997, 2000)           #####
########################################################

##### Load required packages #####
library(vegan)
library(geosphere)
library(mmod)
library(VPLandscapeGenetics)
data(VPLandscapeGenetics)

##### Calculate genetic and geographic distances #####
## Genetic distances
hedrick.SS <- as.matrix(pairwise_Gst_Hedrick(genind.SS))
nei.SS <- as.matrix(pairwise_Gst_Nei(genind.SS))

hedrick.WF <- as.matrix(pairwise_Gst_Hedrick(genind.WF))
nei.WF <- as.matrix(pairwise_Gst_Nei(genind.WF))

## Linearized genetic distances
hedrick.SS.lin <- hedrick.SS / (1-hedrick.SS)
nei.SS.lin <- nei.SS / (1-nei.SS)

hedrick.WF.lin <- hedrick.WF / (1-hedrick.WF)
nei.WF.lin <- nei.WF / (1-nei.WF)

## Geographic distances
geoDist.SS <- as.numeric(distm(latLong.SS, fun=distGeo))
dim(geoDist.SS) <- c(90, 90)
dimnames(geoDist.SS) <- list(1:90, 1:90)

geoDist.WF <- as.numeric(distm(latLong.WF, fun=distGeo))
dim(geoDist.WF) <- c(87, 87)
dimnames(geoDist.WF) <- list(1:87, 1:87)

## Log transform geographic distance
geoDist.SS.log <- log(geoDist.SS)

geoDist.WF.log <- log(geoDist.WF)

##### Perform Mantel tests #####
mantel.hedrick.SS <- mantel(as.dist(hedrick.SS.lin), as.dist(geoDist.SS), permutations = 9999)
mantel.hedrick.log.SS <- mantel(as.dist(hedrick.SS.lin), as.dist(geoDist.SS.log), permutations = 9999)

mantel.nei.SS <- mantel(as.dist(nei.SS.lin), as.dist(geoDist.SS), permutations = 9999)
mantel.nei.log.SS <- mantel(as.dist(nei.SS.lin), as.dist(geoDist.SS.log), permutations = 9999)

mantel.hedrick.WF <- mantel(as.dist(hedrick.WF.lin), as.dist(geoDist.WF), permutations = 9999)
mantel.hedrick.log.WF <- mantel(as.dist(hedrick.WF.lin), as.dist(geoDist.WF.log), permutations = 9999)

mantel.nei.WF <- mantel(as.dist(nei.WF.lin), as.dist(geoDist.WF), permutations = 9999)
mantel.nei.log.WF <- mantel(as.dist(nei.WF.lin), as.dist(geoDist.WF.log), permutations = 9999)
