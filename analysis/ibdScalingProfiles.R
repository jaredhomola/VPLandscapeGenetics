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
library(car)
library(reshape2)
library(VPLandscapeGenetics)
data(VPLandscapeGenetics)

##### Calculate genetic distances #####
nei <- pairwise_Gst_Nei(genind.SS, linearized = TRUE) ## For spotted salamanders
##nei <- pairwise_Gst_Nei(genind.WF, linearized = TRUE) ## For wood frogs

nei.df <- melt(as.matrix(nei), varnames = c("Pop1", "Pop2"))

##### Calculate geographic distances #####
geoDist <- as.dist(distm(latLong.SS, fun=distGeo)) ## For spotted salamanders
##geoDist <- as.dist(distm(latLong.WF, fun=distGeo)) ## For wood frogs

geoDist <- geoDist/1000
attributes(geoDist)$Labels <- attributes(nei)$Labels
geoDist.df <- melt(as.matrix(geoDist), varnames = c("Pop1", "Pop2"))

##### Assemble data frame #####
df <- merge(geoDist.df, nei.df, by = c("Pop1", "Pop2"))
colnames(df) <- c("Pop1", "Pop2", "geoDist", "Gst")

df <- subset(df, geoDist > 0) ## Remove rows with 0 geoDist
df <- df[!duplicated(df$geoDist), ] ## Remove duplicated distances

df <- df[order(df$geoDist), c(1:4)] ## Order data based on geographic distance

########################################################
######             Run analysis loop              ######
########################################################
######    This process should be repeated for     ######
######               each species                 ######
########################################################

i <- 19 # Begins simulation with the first 20 distance points
set.seed(12345) # For reproducibility
Result.CI <- list() # Establish running list of confidence intervals
Result.Beta <- list() # Establish running list of beta values

for(z in 1:nrow(df)) {
  i <- i+1
  x <- df$geoDist[1:i]
  y <- df$Gst[1:i]
  Test <- lm(y ~ x)
  Test.boot <- Boot(Test, R=1000)
  print(paste0("Analysis of first ", i, " distances is complete (", round(max(x), 2), " km)"))
  Result.CI[[length(Result.CI)+1]] = confint(Test.boot, parm=2, level=.95, type="norm")
  Result.Beta[[length(Result.Beta)+1]] = summary(Test.boot, parm=2)
}

## Create dataframe from output
CI.df <- do.call("rbind", lapply(Result.CI, as.data.frame))
Beta.df <- do.call("rbind", lapply(Result.Beta, as.data.frame))
results <- as.data.frame(cbind(Beta.df$original, CI.df$`2.5 %`, CI.df$`97.5 %`, df$geoDist[20:(nrow(Beta.df)+19)]))
colnames(results) <- c("original", "x2.5", "x97.5", "Geo")

###### Plot results ######
# Create appropriate blank plot
plot(1, type="n", xlab="Pairwise distance (km)", ylab="IBD Slope", cex.axis=1.5, cex.lab=1.5, xlim=c(0, 50), ylim=c(-0.002, 0.004))

#UpperCI
UpperCI <- loess(x97.5 ~ Geo, data=results, span=0.05)
UpperCI.smoothed <- predict(UpperCI)

#LowerCI
LowerCI <- loess(x2.5 ~ Geo, data=results, span=0.05)
LowerCI.smoothed <- predict(LowerCI)

#Gray band for CI
polygon(c(rev(results$Geo), results$Geo), c(rev(LowerCI.smoothed), UpperCI.smoothed), col = 'grey85', border = NA)

#Plot beta line
abline(h=0, col="grey55", lwd=2) #Zero line
beta <- loess(original ~ Geo, data=results, span=0.05)
beta.smoothed <- predict(beta)
lines(beta.smoothed, x=results$Geo, col="black", lwd=4)
