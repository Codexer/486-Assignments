# Project 1: Bootstrapping to find the mean
# Confidence interval estimators and comparisons
# Stat 486, Matthew Tantoy (RUID: 163006720, NetID: mat322)

cintBootstrap <- function(dataVec, confidence = 0.95, nboot = 10000) {
	# calculate original statistics and sample size
	nSamp <- length(dataVec)
	origMean <- mean(dataVec)
	origVar <- sum((dataVec - origMean)^2) / (nSamp - 1)

	# generate random samples and save in one giant matrix
	# each row is one resampling
	randInds <- sample.int(nSamp, nboot * nSamp, replace = T)
	sampMat <- matrix(dataVec[randInds], nrow = nboot, ncol = nSamp)

	# studentize samples
	sampMeans <- rowMeans(sampMat)
	sampVars <- rowSums((sampMat - replicate(nSamp, sampMeans))^2) / (nSamp - 1)
	studentizedStats <- (sampMeans - origMean) / sqrt(sampVars / nSamp)

	# with small n, resampling can result in a var of 0
	# which results in Inf and NaN so remove those statistics
	studentizedStats <- studentizedStats[is.finite(studentizedStats)]
	nboot <- length(studentizedStats)

	# use approx to generate "t-score" from studentized stats
	studentizedStats <- sort(studentizedStats)
	interpLimits <- c((1 - confidence) / 2, (1 + confidence) / 2)
	tScores <- approx(x = c(1:nboot) / nboot, y = studentizedStats, xout = interpLimits)$y

	# calculate confidence interval as a vector with two values 
	bounds <- origMean + sqrt(origVar / nSamp) * tScores

	return( bounds )
}

cintCLT <- function(dataVec, confidence = 0.95) {
	cWidth <- qt((1 + confidence) / 2, length(dataVec) - 1) * sqrt(var(dataVec) / length(dataVec))
	return( mean(dataVec) + c(-cWidth, cWidth) )
}

isInInterval <- function(num, interval) {
	return( num < interval[[2]] && num > interval[[1]] )
}

cintTester <- function(sampleSizes = c(3, 10, 30, 100), confidences = c(0.90, 0.95), ntrials = 1000) {
	# NOTE: keeps the same mean and SD for each trial
	testMean <- 0
	testSD <- 1

	for (size in sampleSizes) {
		for (alpha in confidences) {
			countBoot <- 0
			countNorm <- 0
			for (n in 1:ntrials) {
				# generate test set and reference statistic
				# normal distribution
				# testSet <- rnorm(size, testMean, testSD)
				# actualMean <- testMean

				# generate test set and reference statistic
				# lognormal distribution
				testSet <- rlnorm(size, meanlog = testMean, sdlog = testSD)
				actualMean <- exp(testMean + (testSD * testSD) / 2)

				# run intervals with appropriate confidence level
				cintBoot <- cintBootstrap(testSet, confidence = alpha)
				cintNorm <- cintCLT(testSet, confidence = alpha)

				# check coverage rates
				countBoot <- countBoot + isInInterval(actualMean, cintBoot)
				countNorm <- countNorm + isInInterval(actualMean, cintNorm)
			}
			print(sprintf("Tests = Sample size %d, alpha=%f: Bstrap - %d/%d, CLT - %d/%d", 
						size, alpha, countBoot, ntrials, countNorm, ntrials))
		}
	}
}