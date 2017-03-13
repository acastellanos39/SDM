library(raster)
library(fields)
library(RStoolbox)
library(foreach)
PRES <- read.csv("Dipo_pres_1_30.csv", header = T)[, -1]
ABSV <- read.csv("Dipo_background.csv", header = T)[, -1]
MCP <- mcp(SpatialPoints(PRES[, 1:2], proj4string = CRS("+proj=longlat +datum=WGS84")), percent = 100)


###Questions/Thoughts:
##1) If dist() is used on presence points, that may be easier than using the PCA axes, but using up to 95% of PCA axes may be ideal as that would give larger weight to the variables driving the bioclimatic changes within the study extent
##2) How does one determine the scale used for thinning localities without detailed information about the biology/study extent? Do you grab it from distance information about the presence points or from the study extent (leaning more towards study extent)? 
#Perhaps can grab the bioclimatic info for Madagascar where there is a published 10km thinning distance and try to see what distance bioclimatically that represents? 
#^This wasn't really worthwhile, so I didn't include it, but I did write a function to grab the overall mean of the mean distance of a cell and all of its neighbors for any study extent
########################################################

####ENVIRONMENTAL FILTERING BASED ON DISTANCE FUNCTION
#Function that takes presence coordinates, predictor variables, a set distance, and repeats it n times
#Can use the pca.threshold argument to determine at what cumulative proportion of variance to exclude PCA axes from being used to create the distance matrix (I have been using 95%)

balt <- env.filter(PRES[, 1:2], PRED, 0.95, 0.086, 1000) #example

env.filter <- function(presence, predictors, pca.threshold, distance, n) {
	PRES <- extract(predictors, presence)
	PCA <- prcomp(PRES, scale = T)
	INT <- which(summary(PCA)$importance[3, ] > pca.threshold)[1]
	DIST.MAT <- as.matrix(dist(PCA$x[, 1:INT], upper = T))
	REPS <- replicate(n, env.thin(presence, distance), simplify = F)
	list(table(sapply(REPS, length)), sort(table(unlist(REPS)), decreasing = T))	
}

####ENVIRONMENTAL THINNING FUNCTION FOR env.filter()
#Function that is used with the env.filter function to allow for the use of replicate()
#takes the presence coordinates and distance specified from env.filter() and the distance matrix created from  env.filter() to run. Runs a while loop until all acceptable locality records have been assessed for distance

alt <- replicate(1000, env.thin(PRES, 2), simplify = F)

env.thin <- function(presence, distance, distance.matrix = get("DIST.MAT", parent.frame())) {
	KEEP <- seq(1, nrow(presence), 1)
while(length(KEEP) > 1) {
		if(length(KEEP) == nrow(presence)) {
			START <- sample(1:length(KEEP), 1)
			SAMP <- START
			KEEP <- as.numeric(names(which(distance.matrix[, START] > distance)))
			TOSS <- setdiff(seq(1, nrow(presence), 1), KEEP)
		}
		else {
			START <- KEEP[sample(1:length(KEEP), 1)]
			SAMP <- c(SAMP, START)
			KEEP <- as.numeric(names(which(distance.matrix[-TOSS, START] > distance)))
			TOSS <- setdiff(seq(1, nrow(presence), 1), KEEP)
		} 
	}
	SAMP
}


####MEAN DISTANCE FUNCTION
#Function that takes, predictor variables, an extent object (or spatial object that has an extent), a threshold to determine at what point to exclude PCA axes, and a number how many different resolutions you want to investigate
#This function will crop the predictor variables to your region of interest (i.e., to the study extent), run a PCA on each of the raster layers, get rid of extraneous PCA layers based on your given threshold (and determine the weights to be used for the weighted.mean() later), do a sliding window distance calculation for each PCA axis layer, grab the mean for each layer, and then create a weighted overall mean given the proportion of variance explained. 
#Each element of the list is a different resolution. The weighted mean is the first element and the resolution investigated is second. 

alt <- mean.dist(PRED, MCP, 0.95, 5)

mean.dist <- function(predictors, extent, pca.threshold, scale.num) {
	SCALE <- crop(predictors, extent)
	foreach(i = 1:scale.num, .packages = c("foreach", "raster", "RStoolbox")) %do% {
	if(i > 1) {
	SCALE <- aggregate(SCALE, factor = 2)
	}
	PCR <- rasterPCA(SCALE, spca = T)
	INT <- which(cumsum(PCR$model$sdev^2/sum(PCR$model$sdev^2)) > pca.threshold)[1]
	WEIGHT <- (PCR$model$sdev^2/sum(PCR$model$sdev^2))[1:INT]
	MEAN <- foreach(j = 1:INT, .combine = "c") %do% {
		FOCAL <- focal(PCR$map[[j]], w = matrix(1/9, nc = 3, nr = 3), fun = function(x) dist(x))
		cellStats(FOCAL, mean)
		}
	c(weighted.mean(MEAN, WEIGHT), res(SCALE)[1])
	}
}




