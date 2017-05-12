#If one does not know what scale to use, you can use the following function to iterate through a variety of resolutions using the gridSample filtering method. 
#It requires presence coordinates, background coordinates, what k to use in the cross.valid function, predictor variables, how many repetitions to run for each particular scale, what do you want the range of resolutions to be (or more accurately, how many do you want to sample), whether you want to filter in environmental or geographic space, any arguments you so choose for maxent, and what you want inputted in the proj4string argument for spatial filtering.
#The end result will be a list with the same test statistics as from the cross.valid function, but with each column representing a repetition. Each item in the list will be for a different resolution.

range.calc <- function(presence, background, k, predictors, n, range, type = c("environmental", "spatial"), args = NULL, coord.ref = "+proj=longlat +datum=WGS84") {
	if(type = "spatial") {
	PRES <- presence
	RAST <- raster(SpatialPoints(presence, proj4string = CRS(coord.ref)))
	foreach(i = 1:range, .packages = c("dismo", "foreach")) %dopar% {
		res(RAST) <- res(predictors)*i
		SELEC <- presence
		if (i > 1) {
		RET <- gridSample(PRES, RAST, n = 1)
		SELEC <- PRES[rownames(RET), ]
		}
		STAT <- replicate(n, cross.valid(SELEC, background, k, predictors), simplify = F)
		STATS <- list(c(dim(SELEC)[1], res(RAST)[1]), sapply(STAT, '[[', 1), sapply(STAT, '[[', 2))
		STATS
	}	
	}
	if(type = "environmental") {
	PRES <- cbind(presence, extract(predictors, presence))
	PCA <- prcomp(PRES[, -1:-2], scale = T)
	PRES <- cbind(PRES, PCA$x[, 1:2])
	RAST <- raster(extent(range(PRES$PC1), range(PRES$PC2)) + 0.1)
	foreach(i = 1:range, .packages = c("dismo", "foreach")) %dopar% {
		res(RAST) <- res(predictors)*i
		SELEC <- presence
		if(i > 1) {
		RET <- gridSample(PRES[, c("PC1", "PC2")], RAST, n = 1)
		SELEC <- PRES[rownames(RET), 1:2]	
		}
		STAT <- replicate(n, cross.valid(SELEC, background, k, predictors), simplify = F)
		STATS <- list(c(dim(SELEC)[1], res(RAST)[1]), sapply(STAT, '[[', 1), sapply(STAT, '[[', 2))
		STATS
		}
	}
}

STATS <- range.calc(PRES[, 1:2], ABSV[, 1:2], 5, PRED, 100, 10, type = "environmental")


