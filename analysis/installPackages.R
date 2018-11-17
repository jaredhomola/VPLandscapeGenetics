#########################################################
#####         Install packages needed for           #####
#####         VPLandscapeGenetics analyses          #####
#########################################################

list.of.packages <- c("reshape2", "mmod",
                      "geosphere", "glmulti",
                      "hierfstat", "adegenet",
                      "quantreg", "vegan",
                      "tidyverse", "car")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)


