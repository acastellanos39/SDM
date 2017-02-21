#takes presence coordinates, predictor variables, and a scale to filter at and runs the gridSample function at said scale "iter" number of times for the first two PC axes. Grabs the number of samples found highest in all of the runs.
#TO DO List
##check if we even need the iterations component of this function
###does it make a noticeable difference in the model depending on which presence points are chosen (I guess it would depend on the scale, larger scale would allow more difference in the filtered points). In which case, would the gridSample method be best, or would a distance based measure (e.g., something akin to thin from the spThin package) be more appropriate?
##if necessary, write something more to deal with ties (two points chosen half the time in all iterations), tentatively using a distance measure (see last rows 22 and 23 for my start before I balked at going down this route too much)
env.filter <- function(presence.coords, predictors, iter, scale) {
PRES <- cbind(presence.coords, extract(predictors, presence.coords))
PCA <- prcomp(PRES[, -1:-2], scale = T)
RAS <- raster(extent(range(PCA$x[, 1]), range(PCA$x[, 2])) + 0.1)
res(RAS) <- scale
NUM <- replicate(iter, gridSample(as.data.frame(PCA$x[, 1:2]), RAS, n = 1), simplify = F)
NSAMP <- table(sapply(NUM, dim)[1, ])
TABLE <- table(unlist(lapply(NUM, rownames)))
SELEC <- as.numeric(names(sort(TABLE, decreasing = T)[1:names(NSAMP)]))
print(paste("Number of filtered samples:", names(NSAMP)))
presence.coords[SELEC, ]
}

SCALE <- seq(30, 3000, 30)/3600 #can adjust to test or whatever
#alt <- sapply(SCALE, function(x) env.filter(PRES[, 1:2], PRED, 100, x), simplify = F)

#ROWN <- as.numeric(names(which(TABLE == 100/2)))
#ROWM <- sapply(ROWN, function(X) mean(dist(rbind(PCA$x[X, 1:2], PCA$x[, 1:2]))))


#using the env.filter function, can iterate through each chosen resolution and compare based on number of samples kept, AUC, and omission rate. Takes presence coordinates, background coordinates, predictor variables, and a vector of chosen scales to compare.
#TO DO list
##possibly alter the cross validation code here into a function (I already mostly have it done somewhere else, but may edit it to give more information) so I could use a bunch of apply functions to perhaps speed things up(?)
###the maxent function can technically do 5-fold cross validation on its own through an argument, but I didn't trust it because it was giving me weird results with the number of background points entered or something similar
alt <- find.scale(PRES[, 1:2], ABSV[, 1:2], PRED, SCALE[c(5, 10)])
 
find.scale <- function(presence.coords, background, predictors, scales, iter = 1000, k = 5) {
	MCP <- mcp(SpatialPoints(presence.coords[, 1:2], proj4string = CRS("+proj=longlat +datum=WGS84")), percent = 100)
	EXT <- extent(c(MCP@bbox[, 1] - 1, MCP@bbox[, 2] + 1)[c(1, 3, 2, 4)])
	groupb <- kfold(background, k)
	for(i in 1:length(scales)) {
		env.coords <- env.filter(presence.coords, predictors, iter, scale = scales[i])
		group <- kfold(env.coords, k)
		for(j in 1:k) {
		PTRAIN <- env.coords[group != j, ] 
		PTEST <- env.coords[group == j, ] 
		BTRAIN <-background[groupb != j, ] 
		BTEST <- background[groupb == i, ]
		m <- maxent(predictors, PTRAIN, BTRAIN)
		me <- evaluate(PTEST, BTEST, m, predictors)
		if (j == 1) {
		stat <- c(me@auc, abs(me@auc - 	m@results["Training.AUC", ]), abs(sum(me@presence > as.data.frame(m@results)["Minimum.training.presence.logistic.threshold", ])/length(me@presence) - 1))
}
		else {
		stat <- rbind(stat, c(me@auc, abs(me@auc - m@results["Training.AUC", ]), abs(sum(me@presence > as.data.frame(m@results)["Minimum.training.presence.logistic.threshold", ])/length(me@presence) - 1)))
}

		}
		if(i == 1) {
			stats <- as.vector(c(scales[i], dim(env.coords)[1], apply(stat, 2, mean)))
		}
		else {
			stats <- rbind(stats, as.vector(c(scales[i], dim(env.coords)[1], apply(stat, 2, mean))))
		}
	}
	colnames(stats) <- c("Scale", "Num.Coords", "AUC.Test", "AUC.diff", "OmissionRate")
	stats
}


