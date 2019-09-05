########################################################
######         Create Mantel correlograms         ######
########################################################

##### Load required packages #####
library(vegan)
library(mmod)
library(VPLandscapeGenetics)
library(tidyverse)
library(geosphere)
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

###### Run test and produce default plot ######
break.pts <- c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360)
##SS
cor.SS <- mantel.correlog(nei.SS, D.geo=geoDist.SS, mult = "BH", nperm = 1000, break.pts = break.pts, cutoff = TRUE)
plot(cor.SS, alpha=0.05)
SS.results <- as.data.frame(cor.SS$mantel.res)

##WF
cor.WF <- mantel.correlog(nei.WF, D.geo=geoDist.WF, mult = "BH", nperm = 1000, break.pts = break.pts, cutoff = TRUE)
plot(cor.WF)
WF.results <- as.data.frame(cor.WF$mantel.res)

##### Publication plot #####
## Assign pch based on corrected p value using a new column
SS.results$pch <- 1
colnames(SS.results)[5] <- "p.corrected"
SS.results$pch[SS.results$p.corrected <= 0.05] <- 16

WF.results$pch <- 2
colnames(WF.results)[5] <- "p.corrected"
WF.results$pch[WF.results$p.corrected <= 0.05] <- 17

## Generate plot
plot(SS.results$class.index, SS.results$Mantel.cor, ylim=c(-0.08,0.20), xlim=c(20,170), pch=SS.results$pch, font.lab = 2,
      xlab = "Distance class (km)", ylab = "Mantel correlation coefficient", cex=1.5, cex.lab = 1.4)
points(WF.results$class.index, WF.results$Mantel.cor, pch=WF.results$pch, cex=1.5)
lines(WF.results$class.index, WF.results$Mantel.cor)
lines(SS.results$class.index, SS.results$Mantel.cor)
legend(130, 0.225, c("Wood frogs", "Spotted salamanders"),
       pch = c(17, 16), box.lty = 0, bg="transparent", cex = 1.2, pt.cex = 2, y.intersp = .5)

