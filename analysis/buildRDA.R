library(adegenet)

getwd()
envData.SS <- read.delim("G:/My Drive/Homola et al. - Amphibian landscape genetics/R_code/VPLandscapeGenetics/extData/envData.SS.txt", sep="\t", header = TRUE)
latLong.SS <- read.delim("G:/My Drive/Homola et al. - Amphibian landscape genetics/R_code/VPLandscapeGenetics/extData/SS_lat_long.txt", sep="\t", header = FALSE)
genind.SS <- read.genepop("G:/My Drive/Homola et al. - Amphibian landscape genetics/R_code/VPLandscapeGenetics/extData/SS - Filtered.gen", ncode = 3)
envData.SS <- envData.SS[1:10]
latLong.SS <- latLong.SS[2:3]


envData.WF <- read.delim("G:/My Drive/Homola et al. - Amphibian landscape genetics/R_code/VPLandscapeGenetics/extData/envData.WF.txt", sep="\t", header = TRUE)
latLong.WF <- read.delim("G:/My Drive/Homola et al. - Amphibian landscape genetics/R_code/VPLandscapeGenetics/extData/WF_lat_long.txt", sep="\t", header = FALSE)
genind.WF <- read.genepop("G:/My Drive/Homola et al. - Amphibian landscape genetics/R_code/VPLandscapeGenetics/extData/WF - Filtered.gen", ncode = 3)
latLong.WF <- latLong.WF[2:3]

save(envData.SS, latLong.SS, genind.SS, envData.WF, latLong.WF, genind.WF, file = "G:/My Drive/Homola et al. - Amphibian landscape genetics/R_code/VPLandscapeGenetics/data/WFUrbanAdaptation.rda")
