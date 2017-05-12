#performs k-fold cross validation when given presence coordinates, background coordinates, k, predictor variables, and any maxent arguments (if you so choose). 
#This function will grab the mean AUC, difference in training and test AUC, omission rate (using the minimum presence threshold and 10% presence threshold), and true skills statistic (using the same thresholds). Likewise, it grabs the mean percent contribution and permutation importance of each predictor variable. The mean is calculated from each of the k runs. 

cross.valid <- function(presence, background, k, predictors, args = NULL) {
	#MCP <- mcp(SpatialPoints(presence, proj4string = CRS("+proj=longlat +datum=WGS84")), percent = 100)
	#EXT <- extent(c(MCP@bbox[, 1] - 1, MCP@bbox[, 2] + 1)[c(1, 3, 2, 4)])
	group <- list(kfold(presence, k), kfold(background, k))
	DATA <- lapply(seq(1, k, 1), function (x) list(presence[which(group[[1]] == x), ], presence[which(group[[1]] !=  x), ], background[which(group[[2]] == x), ], background[which(group[[2]] != x), ]))
	STAT <- foreach(i = 1:k, .packages= 'dismo') %dopar% {
	m <- maxent(predictors, DATA[[i]][[1]], DATA[[i]][[3]], args = args)
	me <- evaluate(DATA[[i]][[2]], DATA[[i]][[4]], m, predictors)
	#mp <- predict(predictors, m, ext = EXT, progress = '')
	stat <- c(me@auc, abs(me@auc - m@results["Training.AUC", ]), abs(sum(me@presence > as.data.frame(m@results)["Minimum.training.presence.logistic.threshold", ])/length(me@presence) - 1), abs(sum(me@presence > as.data.frame(m@results)["X10.percentile.training.presence.logistic.threshold", ])/length(me@presence) - 1), (me@TPR[which(me@t > as.data.frame(m@results)["Minimum.training.presence.logistic.threshold", ])[1]] + me@TNR[which(me@t > as.data.frame(m@results)["Minimum.training.presence.logistic.threshold", ])[1]] - 1), (me@TPR[which(me@t > as.data.frame(m@results)["X10.percentile.training.presence.logistic.threshold", ])[1]] + me@TNR[which(me@t > as.data.frame(m@results)["X10.percentile.training.presence.logistic.threshold", ])[1]] - 1))
	pred.stat <- c(m@results[paste(names(predictors), ".contribution", sep = ""), ], m@results[paste(names(predictors), ".permutation.importance", sep = ""), ])
	list(stat, pred.stat)
}
	STATS <- apply(t(sapply(STAT, '[[', 1)), 2, mean)
	names(STATS) <- c("AUC", "AUC.diff", "ORate.MPT", "ORate.10PT", "TSS.MPT", "TSS.10PT")
	PRED.STAT <- apply(sapply(STAT, '[[', 2), 1, mean)
	list(STATS, PRED.STAT)
}

CROSS <- cross.valid(PRES[, 1:2], ABSV[, 1:2], 5, PRED)
