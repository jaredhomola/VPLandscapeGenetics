###########################################################
#####     Estimation of IBD scaling profiles to      ######
#####  investigate relationships between regression  ######
#####   beta values and distance class for pairwise  ######
#####  measures of population differentiation and    ######
#####               geographic distances             ######
###########################################################
#####     Script analyzes one species at a time.     ######
#####   Repeat it using the proper input files for   ######
#####                each species                    ######
###########################################################

##### Load required packages #####
library(mmod)
library(geosphere)
library(VPLandscapeGenetics)
library(reshape2)
library(quantreg)
data(VPLandscapeGenetics)

##### Calculate genetic distances #####
## Spotted salamanders
nei.SS <- pairwise_Gst_Nei(genind.SS, linearized = TRUE)
hedrick.SS <- pairwise_Gst_Hedrick(genind.SS, linearized = TRUE)
nei.SS.df <- melt(as.matrix(nei.SS), varnames = c("Pop1", "Pop2"))
hedrick.SS.df <- melt(as.matrix(hedrick.SS), varnames = c("Pop1", "Pop2"))

## Wood frogs
nei.WF <- pairwise_Gst_Nei(genind.WF, linearized = TRUE)
hedrick.WF <- pairwise_Gst_Hedrick(genind.WF, linearized = TRUE)
nei.WF.df <- melt(as.matrix(nei.WF), varnames = c("Pop1", "Pop2"))
hedrick.WF.df <- melt(as.matrix(hedrick.WF), varnames = c("Pop1", "Pop2"))

##### Calculate geographic distances #####
## Spotted salamanders
geoDist.SS <- as.dist(distm(latLong.SS, fun=distGeo))
geoDist.SS <- geoDist.SS/1000
attributes(geoDist.SS)$Labels <- attributes(nei.SS)$Labels
geoDist.SS.df <- melt(as.matrix(geoDist.SS), varnames = c("Pop1", "Pop2"))

## Wood frogs
geoDist.WF <- as.dist(distm(latLong.WF, fun=distGeo))
geoDist.WF <- geoDist.WF/1000
attributes(geoDist.WF)$Labels <- attributes(nei.WF)$Labels
geoDist.WF.df <- melt(as.matrix(geoDist.WF), varnames = c("Pop1", "Pop2"))


##### Assemble data frames for analysis #####
## Spotted salamanders
SS.df <- merge(geoDist.SS.df, nei.SS.df, by = c("Pop1", "Pop2"))
SS.df <- merge(SS.df, hedrick.SS.df, by = c("Pop1", "Pop2"))
colnames(SS.df) <- c("Pop1", "Pop2", "geoDist", "Nei", "Hedrick")
SS.df <- subset(SS.df, geoDist > 0) ## Remove rows with 0 geoDist
SS.df <- SS.df[!duplicated(SS.df$geoDist), ] ## Remove duplicated geographic distances

WF.df <- merge(geoDist.WF.df, nei.WF.df, by = c("Pop1", "Pop2"))
WF.df <- merge(WF.df, hedrick.WF.df, by = c("Pop1", "Pop2"))
colnames(WF.df) <- c("Pop1", "Pop2", "geoDist", "Nei", "Hedrick")
WF.df$geoDist2 <- WF.df$geoDist^2
WF.df <- subset(WF.df, geoDist > 0) ## Remove rows with 0 geoDist
WF.df <- WF.df[!duplicated(WF.df$geoDist), ] ## Remove duplicated geographic distances



###### Set up equations ######
## Spotted salamanders
nei.SS.mod <- lm(Nei ~ geoDist, data = SS.df)
nei.SS.log.mod <- lm(Nei ~ log(geoDist), data = SS.df)
hedrick.SS.mod <- lm(Hedrick ~ geoDist, data = SS.df)
hedrick.SS.log.mod <- lm(Hedrick ~ log(geoDist), data = SS.df)

## Wood frogs
nei.WF.mod <- lm(Nei ~ geoDist, data = WF.df)
nei.WF.log.mod <- lm(Nei ~ log(geoDist), data = WF.df)
hedrick.WF.mod <- lm(Hedrick ~ geoDist, data = WF.df)
hedrick.WF.log.mod <- lm(Hedrick ~ log(geoDist), data = WF.df)

WF.Nei <- WF.df$Nei
WF.geoDist <- WF.df$geoDist
WF.geoDist2 <- WF.df$geoDist2
nei.WF.quadMod <- lm(WF.Nei ~ WF.geoDist + WF.geoDist2) ## Wood frog quadratic model

###### Set up quantile equations ######
## Spotted salamanders
quant95.SS.mod <- rq(Nei ~ geoDist, data=SS.df, tau=.95)
quant05.SS.mod <- rq(Nei ~ geoDist, data=SS.df, tau=.05)

## Wood frogs
quant95.WF.mod <- rq(Nei ~ geoDist + geoDist2, data=WF.df, tau=.95)
quant05.WF.mod <- rq(Nei ~ geoDist + geoDist2, data=WF.df, tau=.05)

###### Perform plotting ######
###### Note: Only Gst plotted due to high similarity between Gst and G''st Mantel tests ######
### Spotted salamanders
plot(SS.df$geoDist, SS.df$Nei, xlab="Pairwise distance (km)", ylab=expression('G'[ST]*'/(1-G'[ST]*')'), pch=20, col="grey50", cex.axis=1.75, cex.lab=1.75, xlim=c(-0.10, 350))
abline(nei.SS.mod, lwd=3, col="black")
abline(quant95.SS.mod, lwd=3, col="black", lty=2)
abline(quant05.SS.mod, lwd=3, col="black", lty=2)

## Wood frogs quadratic model
plot(WF.df$geoDist, WF.df$Nei, xlab="Pairwise distance (km)", ylab=expression('G'[ST]*'/(1-G'[ST]*')'), pch=20, col="grey50", cex.axis=1.75, cex.lab=1.75, xlim=c(-0.10, 350))
geovalues <- seq(0,350, 1)
predictedcounts <- predict(nei.WF.quadMod, list(WF.geoDist = geovalues, WF.geoDist2 = geovalues^2))
lines(geovalues, predictedcounts, col="black", lwd=3)

predictedcounts.quad95 <- predict(quant95.WF.mod, list(geoDist=geovalues, geoDist2=geovalues^2))
lines(geovalues, predictedcounts.quad95, col="black", lwd=3, lty=2)

predictedcounts.quad05 <- predict(quant05.WF.mod, list(geoDist=geovalues, geoDist2=geovalues^2))
lines(geovalues, predictedcounts.quad05, col="black", lwd=3, lty=2)
