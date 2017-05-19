packages <- c("ape", "geosphere", "sp", "ncf", "mclust", "RStoolbox", "adehabitatHR", "raster", "spdep", "plyr", "data.table", "tidyverse")
#lapply(packages, install.packages)
lapply(packages, library, character.only = T)

setwd("~/Desktop/SDM")

#PRES <- read.csv("Dipo_presence.csv", header = T)[, -1]
PRES <- fread("Aud_presence.csv", header = T)[, -1]
MCP <- mcp(SpatialPoints(PRES[, 1:2], proj4string = CRS("+proj=longlat +datum=WGS84")), percent = 100)
EXT <- extent(c(MCP@bbox[, 1] - 1, MCP@bbox[, 2] + 1)[c(1, 3, 2, 4)])

#Bias
BIAS <- fread("bias_aud_alt.csv", header = T)
#BIAS <- mask[, .(.N), by = .(decimalLongitude, decimalLatitude)]
coordinates(BIAS) <- ~decimalLongitude+decimalLatitude
crs(BIAS) <- crs(MCP)
BIAS <- crop(BIAS, EXT)
BIAS <- as.data.table(BIAS)
DIST <- BIAS[, distm(cbind(decimalLongitude, decimalLatitude))]# <-- obviously a horrible idea

#Distance matrices
DIST <- distm(PRES[, 1:2])#distance matrix from geosphere package; apparently calculates great circle distances
#OR
DIST <- sapply(lapply(seq(1, nrow(PRES), 1), function (x) spDistsN1(as.matrix(PRES[, 1:2]), as.numeric(PRES[x, 1:2]), longlat = T)), rbind)#cobbled together distance matrix using the spDistsN1 function from the sp package to calculate distance in kilometers

#spatial scale
##Main question here is what should z be (to replace rownames(PRES) in Moran.I and DIST in correlog)
MORI <- Moran.I(as.numeric(rownames(PRES)), DIST)
CORL <- correlog(PRES[, 1], PRES[, 2], DIST, increment = 2, latlon = T)
#alt <- correlog(BIAS[1:1000, decimalLongitude], BIAS[1:1000, decimalLatitude], BIAS[1:1000, N], increment = 2, latlon = T)

#environmental space scale
PCA <- prcomp(PRES[, -1:-2], scale = T)
summary(PCA)
THRESH <- which(summary(PCA)$importance[3, ] >= 0.95)[1]
MC <- Mclust(PCA$x[, 1], G = 1:20)
summary(MC)

MGRPS <- sapply(lapply(seq(1, THRESH, 1), function(x) Mclust(PCA$x[, x], G = 1:20)$classification), cbind)#way to weight groupings by PC axis proportion of variance? reasonable or not?
