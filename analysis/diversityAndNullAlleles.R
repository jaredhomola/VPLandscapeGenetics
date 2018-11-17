########################################################
######  Estimate allelic richess, examine Hardy-  ######
######    Weinberg proportions, and detection     ######
######       detection of null alleles            ######
########################################################
######         PopGenReport require Latex         ######
######         (e.g., MiKTeX for Windows)         ######
########################################################

##### Load required packages #####

library(PopGenReport)
library(hierfstat)
library(VPLandscapeGenetics)
data(VPLandscapeGenetics)

##### Calculate PopGenReport's report #####
## Spotted salamander
pgr.SS <- popgenreport(genind.SS, mk.counts = TRUE, mk.hwe = TRUE, mk.null.all = TRUE,
                       mk.allel.rich = TRUE, fname = "PopGenReport_SS",
                       foldername = "pgr_SS", path.pgr = NULL)

## Wood frogs
pgr.wf <- popgenreport(genind.WF, mk.counts = TRUE, mk.hwe = TRUE, mk.null.all = TRUE,
                       mk.allel.rich = TRUE, fname = "PopGenReport_WF",
                       foldername = "pgr_WF", path.pgr = NULL)

##### Calculate hierfstat diversity measures and combine into a dataframe #####
## Spotted salamanders
SS.gd <- basic.stats(genind.SS)
SS.Hs <- as.data.frame(t(SS.gd$Hs))
SS.Hs$means <- rowMeans(SS.Hs)
SS.Fis <- as.data.frame(t(SS.gd$Fis))
SS.Fis$means <- rowMeans(SS.Fis)
SS.ar <- as.data.frame(t(as.data.frame(allelic.richness(genind.SS))))
SS.ar$means <- rowMeans(SS.ar)
SS.gd.all <- as.data.frame(cbind(SS.Hs$means, SS.Fis$means, SS.ar$means[2:91]))
names(SS.gd.all) <- c("Hs", "Fis", "AR")
rownames(SS.gd.all) <- rownames(SS.Hs)

## Wood frogs
WF.gd <- basic.stats(genind.WF)
WF.Hs <- as.data.frame(t(WF.gd$Hs))
WF.Hs$means <- rowMeans(WF.Hs)
WF.Fis <- as.data.frame(t(WF.gd$Fis))
WF.Fis$means <- rowMeans(WF.Fis)
WF.ar <- as.data.frame(t(as.data.frame(allelic.richness(genind.WF))))
WF.ar$means <- rowMeans(WF.ar)
WF.gd.all <- as.data.frame(cbind(WF.Hs$means, WF.Fis$means, WF.ar$means[2:88]))
names(WF.gd.all) <- c("Hs", "Fis", "AR")
rownames(WF.gd.all) <- rownames(WF.Hs)

