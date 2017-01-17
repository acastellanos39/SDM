NC <- nested.cross(2, PRES[, 1:2], ABSV[, 1:2], 3, PREDI, BETA[1:2], FEAT[1:2, ], EXT)

##requires dismo, rJava, ENMeval
##function to perform a repeated, nested K-fold cross validation (inner loop uses a subset of the training data to determine the parameter values used in the outer loop). 
#n == number of repetitions
#presence == coordinates of presence data
#background == coordinates of background data
#k == number of folds to use
#predictors == raster stack of environmental variables
#beta == vector of beta multiplier values (given as "betamultiplier==**" to work within the maxent arguments)
#fclass == a matrix of feature class combinations using all five feature classes
#extent == the geographical extent of the study region
nested.cross <- function(n, presence, background, k, predictors, beta, fclass, extent) {
for(f in 1:n) { #loop for repetitions of cross validation
for(i in 1:k) { #outer loop of nested cross validation
group <- kfold(presence, k) 
groupb <- kfold(background, k)
PTRAIN <- presence[group != i, ]
PTEST <- presence[group == i, ]
BTRAIN <- background[groupb != i, ]
BTEST <- background[groupb == i, ]
groupc <- kfold(PTRAIN, k) 
groupd <- kfold(BTRAIN, k)	
for(j in 1:dim(fclass)[1]) { #loop to run through each feature class
		for(g in 1:length(beta)) { #loop to run through each beta value
			for(h in 1:k) { #inner loop to determine parameter values
			PVALID <- PTRAIN[groupc != h, ] #determine validation subsets
			BVALID <- BTRAIN[groupd != h, ]
			m <- maxent(predictors, PVALID, BVALID, args=c(fclass[j, 1], fclass[j, 2], fclass[j, 3], fclass[j, 4], fclass[j, 5], beta[g]))
			mp <- predict(predictors, m, ext = extent, progress = '')
			if(h == 1) {
			aic <- calc.aicc(get.params(m), presence[, 1:2], mp) #calculates AICc values used to determine parameters for outer loop
			}
		else {
			aic <- rbind(aic, calc.aicc(get.params(m), presence[, 1:2], mp))
		}
			}
			if(j == 1 & g ==1) {
			aic.avg <- apply(aic, 2, mean) #averages the AICc from all kfolds of the inner loop to allow comparison
			aic.sd <- apply(aic, 2, sd)	
			}
			else {
			aic.avg <- rbind(aic.avg, apply(aic, 2, mean))
			aic.sd <- rbind(aic.sd, apply(aic, 2, sd))
			}
		}
	}
fcnames <- unlist(lapply(apply(fclass, 1, grep, pattern = "=true", value = T), function(x) paste(sapply(strsplit(x, ""), '[[', 1), collapse = ""))) #grabs first letter of all feature classes that ==true
aicc <- cbind(as.vector(sapply(fcnames, rep, length(beta))), rep(beta, length(fcnames)), aic.avg, aic.sd) #adds in two columns detailing the beta values and feature class of each run
min.fc <-  round(which.min(aicc[, 3])/length(beta) + (0.51 - 1/length(beta))) #determines the lowest AICc value of the runs
m <- maxent(predictors, PTRAIN, BTRAIN, args = c(fclass[min.fc, 1], fclass[min.fc, 2], fclass[min.fc, 3], fclass[min.fc, 4], fclass[min.fc, 5], as.character(aicc[which.min(aicc[, 3]), 2]))) #runs maxent with these parameter values
me <- evaluate(PTEST, BTEST, m, predictors)
mp <- predict(predictors, m, ext = extent, progress = '')
if (i == 1) {
stats <- c(me@auc, m@results["Training.AUC", ], abs(me@auc - m@results["Training.AUC", ]), abs(sum(me@presence > as.data.frame(m@results)["Minimum.training.presence.logistic.threshold", ])/length(me@presence) - 1), m@results["Minimum.training.presence.logistic.threshold", ], aicc[which.min(aicc[, 3]), 1:2]) #for each outer loop run, gives test AUC, training AUC, the difference between the two, omission rate (using the lowest presence threshold), the lowest presence threshold, and the chosen parameters used
}
else {
stats <- rbind(stats, c(me@auc, m@results["Training.AUC", ], abs(me@auc - m@results["Training.AUC", ]), abs(sum(me@presence > as.data.frame(m@results)["Minimum.training.presence.logistic.threshold", ])/length(me@presence) - 1), m@results["Minimum.training.presence.logistic.threshold", ], aicc[which.min(aicc[, 3]), 1:2]))
}		
}
colnames(stats) <- c("Test.AUC", "Training.AUC", "AUC.diff", "Omission.Rate", "LPT", "FC", "Beta")
if(f == 1) {
	stats.all <- stats
	stats.avg <- c(mean(as.numeric(stats[, 1])), sd(as.numeric(stats[, 1]))) #averages and grabs the standard deviation for the test AUC for the combined outer loop runs
	print(stats.avg) #prints the average and standard deviation
}
else {
	stats.all <- rbind(stats.all, stats)	
	stats.avg <- rbind(stats.avg, c(mean(as.numeric(stats[, 1])), sd(as.numeric(stats[, 1]))))
	print(stats.avg[n, ])
}
write.csv(stats.all, "results_all.csv") #writes all the evaulation statistics for each separate run to a .csv file
write.csv(stats.avg, "results_avg.csv") #writes the average of all runs (at the time of the repetition e.g., for a kfold of 5, 5 for first repetition, 10 for second, 15 for third, etc.) to a .csv file
}
list(stats.all, stats.avg)
}
