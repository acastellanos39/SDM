#This code will run through environmental filters made by every possible combination of predictor variables. It gives a data frame detailing the number of variables used to filter, which variables were used to filter, the number of samples after filtering, and the average evaluation statistics of the five cross validation runs. Before the run, you need to have the PRES, ABSV, and PREDI objects. 

##Notes: If needed, we can grab standard deviation of any of the evaluation statistics as well. I didn't run any of the additional Maxent arguments or test for the best parameters for each run to keep computation time down and because I figure that it is fine to keep everything at default as long as each comparison is at default as well.


#use the stack below to test if needed -- it won't take long to run it with the below variables
#PREDI <- stack(PRED$bio2, PRED$bio4, PRED$bio6, PRED$bio8)

k <- 5
SEL <- list()
SELECT <- list()
groupb <- kfold(ABSV, k)
for(n in 2:length(names(PREDI))) {
	alt <- combn(names(PREDI), n)
	if(nrow(alt) > 2) {
		for(i in 1:dim(alt)[2]) {
			balt <- combn(alt[, i], 2)
			dalt <- PRES
			for(m in 1:dim(balt)[2]) {
				RAST <- raster(extent(range(dalt[balt[1, m]]), range(dalt[balt[2, m]])) + 0.1)
	res(RAST) <- res(PREDI)
	SELE <- gridSample(dalt[balt[, m]], RAST, n = 1)
	calt <- dalt[rownames(SELE), ]
	dalt <- calt
			}
	SELECT[[i]] <- dalt
		}
	}
	else {
	for(i in 1:dim(alt)[2]) {
	RAST <- raster(extent(range(PRES[alt[1, i]]), range(PRES[alt[2, i]])) + 0.1)
	res(RAST) <- res(PREDI)
	SEL[[i]] <- gridSample(PRES[alt[, i]], RAST, n = 1)
	SELECT[[i]] <- PRES[rownames(SEL[[i]]), ]
	}
}
for(j in 1:dim(alt)[2]) {
	colnames(SELECT[[j]])[1:2] <- c("lon", "lat")
	set.seed(7)
	group <- kfold(SELECT[[j]], k)
	for(g in 1:k) {
	PTRAIN <- SELECT[[j]][group != g, 1:2]
	PTEST <- SELECT[[j]][group == g, 1:2]
	BTRAIN <- ABSV[groupb != g, 1:2]
	BTEST <- ABSV[groupb == g, 1:2]
	M <- maxent(PREDI, PTRAIN, BTRAIN) 
	ME <- evaluate(PTEST, BTEST, M, PREDI)
	if (g == 1) {
		stats <- c(ME@auc, M@results["Training.AUC", ], abs(ME@auc - M@results["Training.AUC", ]), abs(sum(ME@presence > as.data.frame(M@results)["Minimum.training.presence.logistic.threshold", ])/length(ME@presence) - 1), M@results["Minimum.training.presence.logistic.threshold", ])
	}
	else {
		stats <- rbind(stats, c(ME@auc, M@results["Training.AUC", ], abs(ME@auc - M@results["Training.AUC", ]), abs(sum(ME@presence > as.data.frame(M@results)["Minimum.training.presence.logistic.threshold", ])/length(ME@presence) - 1), M@results["Minimum.training.presence.logistic.threshold", ]))
	}}
	if(j == 1) {
	SUMTAB <- c(dim(alt)[1], apply(as.matrix(alt[, j]), 2, paste, collapse = " "), dim(SELECT[[j]])[1], apply(stats, 2, mean))	
	}
	else {
		SUMTAB <- rbind(SUMTAB, c(dim(alt)[1], apply(as.matrix(alt[, j]), 2, paste, collapse = " "), dim(SELECT[[j]])[1], apply(stats, 2, mean)))
	}
}
if(n == 2) {
	SUMT <- SUMTAB
}
else {
	SUMT <- rbind(SUMT, SUMTAB)
}
colnames(SUMT) <- c("Num.var", "Var", "Num.samp", "Test.AUC", "Train.AUC", "AUC.diff", "Omission", "LPT")
rownames(SUMT) <- NULL
print(as.data.frame(SUMT))
}


