# nboot=100, ntrials=50 performance notes:
# original var function - ~4.62 sec
# custom myVar2 function (not calculating mean or length) - ~2.52 sec
# custom myMean function (not calculating length) - ~2.12 sec
# change sample to sample.int - ~1.88 sec 

# complete overhaul and vectorize everything - ~0.5 sec

# myVar <- function(dataVec) {
# 	sampMean <- mean(dataVec)
# 	total <- 0
# 	for (num in dataVec) {
# 		total <- total + (num - sampMean) * (num - sampMean)
# 	}

# 	return( total / length(dataVec) )
# }

myVar2 <- function(dataVec, sampMean, sampLen) {
	sum((dataVec - sampMean)^2) / sampLen
}

myMean <- function(dataVec, sampLen) {
	sum(dataVec) / sampLen
}

# mySample <- function(dataVec, sampLen) {
# 	dataVec[floor(runif(sampLen, min = 0, max = sampLen)) + 1]
# }

createProfileData <- function() {
	Rprof()
	cintTester()
	Rprof(NULL)
}

# R doesn't have 64 bit integers (lol), so use a 32-bit xorshift*
# Note: I have no idea what I'm doing with PRNGs
# EDIT: R doesn't accept -2^31 (lol) so shifts may (rarely)
# return NA instead
rState <- as.integer(strtoi('0x2545F491'))
# rMult <- 1332534557L
getRandNum <- function() {
	rState <<- bitwXor(rState, bitwShiftR(rState, 9L))
	rState <<- bitwXor(rState, bitwShiftR(rState, 2L))
	rState <<- bitwXor(rState, bitwShiftL(rState, 7L))
	rState
}

# Remove bias
getRandNumBound <- function(bound) 
{	threshold <- 4294967296 %% bound
	repeat {
		if (getRandNum() >= threshold) {
			return( rState %% bound )
		}
	}
}

cintBootstrap <- function(dataVec, confidence = 0.95, nboot = 10000) {
	nSamp <- length(dataVec)
	origMean <- mean(dataVec)
	origVar <- sum((dataVec - origMean)^2) / (nSamp - 1)

	# randInds <- matrix(sample.int(nSamp, nboot * nSamp, replace = T), nrow = nboot, ncol = nSamp)
	# repInds <- rep(TRUE, nboot)
	# for (i in 2:nSamp) {
	# 	repInds <- repInds & (randInds[,1] == randInds[,i])
	# }
	# while (any(repInds)) {
	# 	randInds[repInds, ] <- sample.int(nSamp, sum(repInds) * nSamp, replace = T)
	# 	repInds <- rep(TRUE, nboot)
	# 	for (i in 2:nSamp) {
	# 		repInds <- repInds & (randInds[,1] == randInds[,i])
	# 	}
	# }
	# sampMat <- matrix(dataVec[randInds], nrow = nboot, ncol = nSamp)

	randInds <- sample.int(nSamp, nboot * nSamp, replace = T)
	sampMat <- matrix(dataVec[randInds], nrow = nboot, ncol = nSamp)

	sampMeans <- rowMeans(sampMat)
	sampVars <- rowSums((sampMat - replicate(nSamp, sampMeans))^2) / (nSamp - 1)
	studentizedStats <- (sampMeans - origMean) / sqrt(sampVars / nSamp)

	# if (any(is.infinite(studentizedStats))) {
	# 	print(sampMat[is.infinite(studentizedStats), ])
	# 	print(sampMeans[is.infinite(studentizedStats)])
	# 	print(sampVars[is.infinite(studentizedStats)])
	# }

	# print(sampMat[1, ])
	# print(head(sampMeans))
	# print(head(sampVars))
	# print(head(studentizedStats))

	# readline()

	studentizedStats <- studentizedStats[is.finite(studentizedStats)]
	nboot <- length(studentizedStats)

	studentizedStats <- sort(studentizedStats)
	interpLimits <- c((1 - confidence) / 2, (1 + confidence) / 2)
	tScores <- approx(x = c(1:nboot) / nboot, y = studentizedStats, xout = interpLimits)$y
	bounds <- origMean + sqrt(origVar / nSamp) * tScores

	# if (any(is.nan(bounds))) {
	# 	print(head(studentizedStats))
	# 	print(interpLimits)
	# 	print(tScores)
	# }

	return( bounds )
}

cintTester <- function(sampleSizes = c(3, 10, 30, 100), confidences = c(0.90, 0.95), ntrials = 1000) {
	testMean <- 0
	testSD <- 1

	for (size in sampleSizes) {
		for (alpha in confidences) {
			countBoot <- 0
			countNorm <- 0
			for (n in 1:ntrials) {
				# generate test set and reference statistic
				# testSet <- rnorm(size, testMean, testSD)
				# actualMean <- testMean

				testSet <- rlnorm(size, meanlog = testMean, sdlog = testSD)
				actualMean <- exp(testMean + (testSD * testSD) / 2)

				cintBoot <- cintBootstrap(testSet, confidence = alpha)
				cintNorm <- cintCLT(testSet, confidence = alpha)

				countBoot <- countBoot + isInInterval(actualMean, cintBoot)
				countNorm <- countNorm + isInInterval(actualMean, cintNorm)

				# print(cintBoot)
				# print(cintNorm)
				# print(actualMean)
				# print(countBoot)
				# print(countNorm)

				# if (is.na(countBoot)) {
				# 	readline()	
				# }
			}
			print(sprintf("Tests = Sample size %d, alpha=%f: Bstrap - %d/%d, CLT - %d/%d", 
						size, alpha, countBoot, ntrials, countNorm, ntrials))
			# readline()
		}
	}
}

timeTest <- function(dataVec) {
	starttime <- Sys.time()
	throwaway <- cintBootstrap(dataVec)
	endtime <- Sys.time()
	return( endtime - starttime )
}