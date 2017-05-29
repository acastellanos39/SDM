packages <- c("tidyverse", "foreach", "raster")
#lapply(packages, install.packages)
lapply(packages, library, character.only = T)

setwd("~/Desktop/SDM")
PRES <- read.csv("Dipo_presence.csv", header = T)[, -1]

###ENVIRONMENTAL FILTERING USING PCA
#the two below functions are similar to what I wrote for EnvFilterDistanceR2.R. The main difference is that the idea behind this pair of functions isn't to use a distance matrix derived from Euclidean distances of PCA axes, but filter similar to the gridSample function (just along one axis at a time rather than two). 
#Again, there are two functions to allow for replicate() to be used. The thin() function will break each axis into a series of bins based on the scale chosen and randomly sample one point from each bin (if possible, since some bins will have no points in them). When moving to the next PC axis for further filtering, only points kept from the last filtering step will be considered. 
thin <- function(pca.data = get("DATA", parent.frame()), scale) {
	DATA <- pca.data
for (i in 1:THRESH) {
	if (i == 1) {
		KEEP <- unlist(map(Filter(length, 	split(rownames(DATA), cut(DATA[, i], seq(min(DATA[, i]), max(DATA[, i]), by = scale)))), sample, 1))
	}
	else {
		DATA <- DATA[KEEP, ]
		KEEP <- unlist(map(Filter(length, split(rownames(DATA), cut(DATA[, i], seq(min(DATA[, i]), max(DATA[, i]), by = scale)))), sample, 1))
	}
}
KEEP
}

#The extraction of values and the PCA is conducted within this function to avoid an unnecessary of waste time and resources. The number of replications (default set here at 1,000) and range of scales to be considered is to designated in this function. The output for this function is a list detailing the scale used, number of samples kept after filtering for each repetition, and the number of times each particular sample was kept after all filtering steps. I used purrr:map() quite a bit because I found it more elegant than lapply in this situation, but I can easily replace it with lapply to reduce the reliance on other packages.
thin.filter <- function(presence, predictors, n = 1000, range) {
	PRES <- presence
	PRESD <- raster::extract(predictors, PRES)
	PCA <- prcomp(PRESD, scale = T)
	THRESH <- which(summary(PCA)$importance[3, ] >= 0.95)[1]
	DATA <- as.data.frame(PCA$x[, 1:THRESH])
	foreach(j = 1:length(range), .packages = "tidyverse") %do% {
	REPS <- replicate(n, thin(DATA, range[j]), simplify = F)
	TABS <- list(range[j], Num.Samples = table(unlist(map(REPS, length))), Filtered.Rownums = sort(table(unlist(map(REPS, as.character))), decreasing = T))
	TABS	
	}	
}

thin.filter(PRES[, 1:2], PRED, 10, range = seq(res(PRED)[1], res(PRED)[1]*10, by = res(PRED)[1]))
