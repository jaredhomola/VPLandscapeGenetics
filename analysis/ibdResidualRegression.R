#########################################################
#####   Generate population-wise mean IBD residual  #####
##### values and analyze environmental associations #####
#########################################################

##### Load required packages #####
library(reshape2)
library(mmod)
library(geosphere)
library(glmulti)
library(VPLandscapeGenetics)
data(VPLandscapeGenetics)

##### Calculate genetic distances #####
## Spotted salamanders
nei.SS <- pairwise_Gst_Nei(genind.SS, linearized = TRUE)
nei.SS.df <- melt(as.matrix(nei.SS), varnames = c("Pop1", "Pop2"))

## Wood frogs
nei.WF <- pairwise_Gst_Nei(genind.WF, linearized = TRUE)
nei.WF.df <- melt(as.matrix(nei.WF), varnames = c("Pop1", "Pop2"))


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
colnames(SS.df) <- c("Pop1", "Pop2", "geoDist", "Nei")
SS.df <- subset(SS.df, geoDist > 0) ## Remove rows with 0 geoDist
SS.df <- SS.df[!duplicated(SS.df$geoDist), ] ## Remove duplicated geographic distances

WF.df <- merge(geoDist.WF.df, nei.WF.df, by = c("Pop1", "Pop2"))
colnames(WF.df) <- c("Pop1", "Pop2", "geoDist", "Nei")
WF.df$geoDist2 <- WF.df$geoDist^2
WF.df <- subset(WF.df, geoDist > 0) ## Remove rows with 0 geoDist
WF.df <- WF.df[!duplicated(WF.df$geoDist), ] ## Remove duplicated geographic distances


###### Set up equations ######
## Spotted salamanders
nei.SS.mod <- lm(Nei ~ geoDist, data = SS.df)

## Wood frogs (Quadratic model)
nei.WF.mod <- lm(Nei ~ geoDist + geoDist2, data = WF.df) ## Wood frog quadratic model


###### Calculate IBD residuals ######
## Quantify residuals ##
SS.resid.vector <- as.vector(resid(nei.SS.mod))
WF.resid.vector <- as.vector(resid(nei.WF.mod))

## Pair residuals with populations ##
SS.pops.resids1 <- as.data.frame(SS.df[1])
SS.pops.resids2 <- as.data.frame(SS.df[2])

SS.pops.resids1$Resids <- SS.resid.vector
SS.pops.resids2$Resids <- SS.resid.vector
names(SS.pops.resids1)[names(SS.pops.resids1)=="Pop1"] <- "pop"
names(SS.pops.resids2)[names(SS.pops.resids2)=="Pop2"] <- "pop"

SS.pops.resids <- rbind(SS.pops.resids1, SS.pops.resids2)

WF.pops.resids1 <- as.data.frame(WF.df[1])
WF.pops.resids2 <- as.data.frame(WF.df[2])

WF.pops.resids1$Resids <- WF.resid.vector
WF.pops.resids2$Resids <- WF.resid.vector
names(WF.pops.resids1)[names(WF.pops.resids1)=="Pop1"] <- "pop"
names(WF.pops.resids2)[names(WF.pops.resids2)=="Pop2"] <- "pop"

WF.pops.resids <- rbind(WF.pops.resids1, WF.pops.resids2)

## Calculate average residual per population ##
## Mean function
i=1
SSOutput <- list()

for(i in 1:length(unique(SS.pops.resids$pop))){
  iname <- as.character(unique(SS.pops.resids$pop)[i])
  ipop <- subset(SS.pops.resids, SS.pops.resids$pop==iname)
  SSresults <- c(mean(ipop$Resids))
  SSOutput[[length(SSOutput)+1]] = SSresults
  print(iname)
}

as.vector(SSOutput)
SSResult.df <- do.call("rbind", lapply(SSOutput, as.data.frame))
names(SSResult.df)[names(SSResult.df)=="X[[i]]"] <- "Resids"

i=1
WFOutput <- list()

for(i in 1:length(unique(WF.pops.resids$pop))){
  iname <- as.character(unique(WF.pops.resids$pop)[i])
  ipop <- subset(WF.pops.resids, WF.pops.resids$pop==iname)
  WFresults <- c(mean(ipop$Resids))
  WFOutput[[length(WFOutput)+1]] = WFresults
  print(iname)
}

as.vector(WFOutput)
WFResult.df <- do.call("rbind", lapply(WFOutput, as.data.frame))
names(WFResult.df)[names(WFResult.df)=="X[[i]]"] <- "Resids"


###### Analyze relationships between IBD residuals and environmental features ######
## Assemble dataframe
SS.resid.df <- as.data.frame(cbind(SSResult.df$Resids, envData.SS$SS_lightroads,
                                   envData.SS$SS_secondaryroads, envData.SS$SS_primaryroads))
names(SS.resid.df) <- c("meanResid", "lightRoads", "secondaryRoads", "primaryRoads")

WF.resid.df <- as.data.frame(cbind(WFResult.df$Resids, envData.WF$WF_lightroads,
                                   envData.WF$WF_secondaryroads, envData.WF$WF_primaryroads))
names(WF.resid.df) <- c("meanResid", "lightRoads", "secondaryRoads", "primaryRoads")

## Variable selection
vs.SS <- glmulti(meanResid ~ lightRoads + secondaryRoads + primaryRoads, data=SS.resid.df, level=1, fitfunction=glm, crit="aicc")
print(vs.SS)
weightable(vs.SS)

vs.WF <- glmulti(meanResid ~ lightRoads + secondaryRoads + primaryRoads, data=WF.resid.df, level=1, fitfunction=glm, crit="aicc")
print(vs.WF)
weightable(vs.WF)

## Full models
full.SS <- lm(meanResid ~ lightRoads + secondaryRoads + primaryRoads, data=SS.resid.df)
summary(full.SS)

full.WF <- lm(meanResid ~ lightRoads + secondaryRoads + primaryRoads, data=WF.resid.df)
summary(full.WF)

## Univariate models
uni.SS <- lm(meanResid ~ lightRoads, data=SS.resid.df)
summary(uni.SS)

uni.WF <- lm(meanResid ~ lightRoads, data=WF.resid.df)
summary(uni.WF)

## Plot univariate models
plot(SS.resid.df$lightRoads, SS.resid.df$meanResid, xlab="Length of light roads within 1km (m)", ylab="Mean population IBD residual", main="Spotted salamanders", pch=20, cex = 2, col="grey50", cex.axis=1.75, cex.lab=1.75, ylim=c(-0.01, 0.024), xlim=c(0, 16000))
abline((lm(SS.resid.df$meanResid ~ SS.resid.df$lightRoads)), lwd=3, col="black")

plot(WF.resid.df$lightRoads, WF.resid.df$meanResid, xlab="Length of light roads within 1km (m)", ylab="Mean population IBD residual", main="Spotted salamanders", pch=20, cex = 2, col="grey50", cex.axis=1.75, cex.lab=1.75, ylim=c(-0.01, 0.024), xlim=c(0, 16000))
abline((lm(WF.resid.df$meanResid ~ WF.resid.df$lightRoads)), lwd=3, col="black")
