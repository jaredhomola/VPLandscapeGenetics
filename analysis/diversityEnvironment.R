#########################################################
#####  Analysis of environmental correlations with  #####
#####          measures of genetic diversity        #####
#########################################################

##### Load required packages #####
library(hierfstat)
library(VPLandscapeGenetics)
data(VPLandscapeGenetics)

##### Calculate diversity measures and combine into a dataframe
## Spotted salamanders
SS.gd <- basic.stats(genind.SS)
SS.Hs <- as.data.frame(t(SS.gd$Hs))
SS.Hs$means <- rowMeans(SS.Hs)
SS.ar <- as.data.frame(t(as.data.frame(allelic.richness(genind.SS))))
SS.ar$means <- rowMeans(SS.ar)
SS.gd.all <- as.data.frame(cbind(SS.Hs$means, SS.ar$means[2:91]))
names(SS.gd.all) <- c("Hs", "AR")
rownames(SS.gd.all) <- rownames(SS.Hs)

## Wood frogs
WF.gd <- basic.stats(genind.WF)
WF.Hs <- as.data.frame(t(WF.gd$Hs))
WF.Hs$means <- rowMeans(WF.Hs)
WF.ar <- as.data.frame(t(as.data.frame(allelic.richness(genind.WF))))
WF.ar$means <- rowMeans(WF.ar)
WF.gd.all <- as.data.frame(cbind(WF.Hs$means, WF.ar$means[2:88]))
names(WF.gd.all) <- c("Hs", "AR")
rownames(WF.gd.all) <- rownames(WF.Hs)

##### Environment-Genetic diversity regression models #####
SS.Hs.mod <- lm(SS.gd.all$Hs ~ envData.SS$SS_lightroads + envData.SS$SS_secondaryroads + envData.SS$SS_primaryroads)
WF.Hs.mod <- lm(WF.gd.all$Hs ~ envData.WF$WF_lightroads + envData.WF$WF_secondaryroads + envData.WF$WF_primaryroads)

SS.AR.mod <- lm(SS.gd.all$AR ~ envData.SS$SS_lightroads + envData.SS$SS_secondaryroads + envData.SS$SS_primaryroads)
WF.AR.mod <- lm(WF.gd.all$AR ~ envData.WF$WF_lightroads + envData.WF$WF_secondaryroads + envData.WF$WF_primaryroads)

##### Isolation-Genetic diversity regression models #####
##### Note: Isolation metrics saved in envData files #####
##### but calculated using ibdResidualRegress.R script #####

SS.HsIso.mod <- lm(SS.gd.all$Hs ~ envData.SS$SS_resids)
WF.HsIso.mod <- lm(WF.gd.all$Hs ~ envData.WF$WF_resids)

SS.ARIso.mod <- lm(SS.gd.all$AR ~ envData.SS$SS_resids)
WF.ARIso.mod <- lm(WF.gd.all$AR ~ envData.WF$WF_resids)

##### Plot isolation-Genetic diversity regression models #####

plot(envData.SS$SS_resids, SS.gd.all$Hs, xlab="Population average IBD residual", ylab="Expected heterozygosity", main="Spotted salamanders", pch=20, cex.axis=1.5, cex.lab=1.5)
abline(SS.HsIso.mod, lwd=3)

plot(envData.WF$WF_resids, WF.gd.all$Hs, xlab="Population average IBD residual", ylab="Expected heterozygosity", main="Wood frogs", pch=20, cex.axis=1.5, cex.lab=1.5)
abline(WF.HsIso.mod, lwd=3)


plot(envData.SS$SS_resids, SS.gd.all$AR, xlab="Population average IBD residual", ylab="Allelic richness", main="Spotted salamanders", pch=20, cex.axis=1.5, cex.lab=1.5)
abline(SS.ARIso.mod, lwd=3)

plot(envData.WF$WF_resids, WF.gd.all$AR, xlab="Population average IBD residual", ylab="Allelic richness", main="Wood frogs", pch=20, cex.axis=1.5, cex.lab=1.5)
abline(WF.ARIso.mod, lwd=3)


