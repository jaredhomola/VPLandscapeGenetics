########################################################
####   Distance-based Redundancy Analyses (dbRDA)   ####
########################################################

##### Load required packages #####
library(dplyr)
library(vegan)
library(mmod)
library(VPLandscapeGenetics)
data(VPLandscapeGenetics)

##### Calculate genetic distances #####
nei.SS <- as.matrix(pairwise_Gst_Nei(genind.SS))
nei.WF <- as.matrix(pairwise_Gst_Nei(genind.WF))

#####  Perform dbRDA  ######
## Evaluate full models
SS.dbRDA <- capscale(nei.SS ~ envData.SS$SS_lightroads + envData.SS$SS_secondaryroads + envData.SS$SS_primaryroads + Condition(latLong.SS[,1] + latLong.SS[,2]))
SS.dbRDAr <- anova.cca(SS.dbRDA, permutations = 10000)
anova(SS.dbRDA, by="term", permu=10000) ## test for sig. environ. variables

WF.dbRDA <- capscale(nei.WF ~ envData.WF$WF_lightroads + envData.WF$WF_secondaryroads + envData.WF$WF_primaryroads + Condition(latLong.WF[,1] + latLong.WF[,2]))
WF.dbRDAr <- anova.cca(WF.dbRDA, permutations = 10000)
anova(WF.dbRDA, by="term", permu=10000) ## test for sig. environ. variables

## Backward optimized models (only light roads retained)
SS.dbRDA2 <- capscale(nei.SS ~ envData.SS$SS_lightroads + Condition(latLong.SS[,1] + latLong.SS[,2]))
SS.dbRDAr2 <- anova.cca(SS.dbRDA2, permutations = 10000)

WF.dbRDA2 <- capscale(nei.WF ~ envData.WF$WF_lightroads + Condition(latLong.WF[,1] + latLong.WF[,2]))
WF.dbRDAr2 <- anova.cca(WF.dbRDA2, permutations = 10000)

## Variance partitioning
mod.SS2 <- varpart(nei.SS, ~envData.SS$SS_lightroads, ~latLong.SS[,1] + latLong.SS[,2])
plot(mod.SS2)
mod.WF2 <- varpart(nei.WF, ~envData.WF$WF_lightroads, ~latLong.WF[,1] + latLong.WF[,2])
plot(mod.WF2)
