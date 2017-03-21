regpluspress <- function(x, y) {
	ls.str <- lsfit(x, y)
	ls.str$leverage <- hat(x)
	ls.str$press <- sum( (ls.str$resid / (1 - ls.str$leverage))^2 )
	return( ls.str )
}

matrix.2ndorder.make <- function(x, only.quad=F) {
	x0<-x
	dimn<-dimnames(x)[[2]] #extract the names of the variables
	num.col<-length(x[1,]) # how many columns

	for(i in 1:num.col) {
		# if we are doing all 2nd order
		if(!only.quad) {
			for(j in i:num.col) {
				x0<-cbind(x0,x[,i]*x[,j])
				dimn<-c(dimn,paste(dimn[i],dimn[j],sep="*"))#create interaction dimnames
			}
		} else {
			#in here only if doing only squared terms
			x0<-cbind(x0,x[,i]*x[,i])
			dimn<-c(dimn,paste(dimn[i],"2",sep="")) # squared dimension names
		}
	}

	dimnames(x0)[[2]]<-dimn #assign the new dimension names
	x0
}

pressReg <- function(xmat0, yvec, print.ls = FALSE, debug.plot = FALSE, quadChoice = FALSE) {
	# add quadratic variables
	xmat <- matrix.2ndorder.make(xmat0, only.quad = quadChoice)

	# run leaps
	leaps.str <- leaps(xmat, yvec)

	# rank data by Cp (calculate Cp - p and order)
	z1 <- leaps.str$Cp - leaps.str$size
	o1 <- order(z1)
	parvec <- (leaps.str$which[o1, ])[1, ]
	xsubset <- xmat[, parvec]
	Cpval <- z1[o1][[1]]

	# perform regression
	ls.str0 <- regpluspress(xsubset, yvec)

	# use regression to make prediction
	coef <- ls.str0$coef
	ypred <- xsubset %*% c(coef[-1]) + coef[1]

	# calculate MPSE (debug only)
	npar <- sum(parvec)
	mpse <- ls.str0$press / (length(yvec) - (npar + 1))
	MPSE <- floor(1000 * mpse / 1000)

	if (print.ls) {
		ls.print(ls.str0)

		print(head(ypred))
		print(sprintf('Press = %f', ls.str0$press))
		print(sprintf('MPSE = %f', mpse))
		print(sprintf('Cp - p = %f', Cpval))
	}

	if (debug.plot) {
		# setup plotting area
		par(mfrow = c(1, 2))

		plot(ypred, yvec, 
			main = sprintf("MPSE = %f, Cp - p = %f", 
							MPSE, Cpval)
		)
		plot(ypred, ls.str0$resid)
	}

	# return variables selected and regression coefficients
	return( list(xselect = parvec, ls.str = ls.str0) )
}

predictAtPoint <- function(xpred0, xmat0, y, quadChoice = FALSE) {
	# variable selection and regression
	regResult <- pressReg(xmat0, y, quadChoice = quadChoice)
	ls.str <- regResult$ls.str

	# perform quadratic expansion and variable selection on prediction row
	xpredFull <- matrix.2ndorder.make(xpred0, only.quad = quadChoice)
	xpred <- xpredFull[, regResult$xselect]

	# calculate prediction
	ypred <- ls.str$coef %*% c(1,xpred)

	# use ls.diag to extract covariance matrix
	ycov<-ls.diag(ls.str)$cov.unscaled
	# use ls.diag to extract std deviation
	std.dev<-ls.diag(ls.str)$std.dev

	# calculate degrees of freedom using ycov
	df <- length(y) - length(diag(ycov))

	# variance of data around line
	v1<-std.dev^2
	# variance of prediction
	vpred<-v1*c(1,xpred) %*% ycov %*% c(1,xpred)
	
	#output prediction and standard deviation of prediction
	list(pred=ypred, sd=sqrt(vpred), df=df)
}

pintPerfect <- function(xpred0, xmat0, y, quadChoice = TRUE, confidence = 0.95) {
	predResult <- predictAtPoint(xpred0, xmat0, y, quadChoice = quadChoice)
	width <- qt((1 + confidence) / 2, predResult$df) * predResult$sd
	return( predResult$pred + c(-width, width) )
}

pintBootstrap <- function(xpred0, xmat0, y, nboot = 10000, confidence = 0.95, quadChoice = TRUE) {
	# perform normal prediction calculations
	origResult <- predictAtPoint(xpred0, xmat0, y, quadChoice = quadChoice)
	origPred <- origResult$pred
	origSd <- origResult$sd

	# initialize memory, get number of samples
	nSamp <- length(y)
	studStats <- c(1:nboot)

	# bootstrap loop
	for (i in 1:nboot) {
		# print("bootstrapping at")
		# print(i)

		# resample random cars from the dataset with replacement
		randInds <- sample.int(nSamp, nSamp, replace = T)
		currX <- xmat0[randInds, , drop = FALSE]
		currY <- y[randInds]

		# calculate t-score from that random selection
		currResult <- predictAtPoint(xpred0, currX, currY, quadChoice = quadChoice)
		studStats[i] <- (currResult$pred - origPred) / (currResult$sd)  

		# readline()
	}

	# if somehow sd = 0,
	# some stats are Inf / NaN so remove those statistics
	# (probably not necessary, but has a negligible effect on runtime)
	studStats <- studStats[is.finite(studStats)]
	nboot <- length(studStats)

	# use approx to generate "t-score" from studentized stats
	studStats <- sort(studStats)
	interpLimits <- c((1 - confidence) / 2, (1 + confidence) / 2)
	tScores <- approx(x = c(1:nboot) / nboot, y = studStats, xout = interpLimits)$y

	# print(tScores)
	# readline()

	# calculate confidence interval as a vector with two values 
	bounds <- origPred + origSd * tScores

	return( bounds )
}

main <- function() {
	# load data set
	library(ISLR)
	# load leaps calculation
	library(leaps)

	# note: removes origin variable
	autoMat <- as.matrix(Auto[, 1:7])
	x.auto <- autoMat[, -1]
	y.auto <- autoMat[, 1]

	# take first row, cut weight in half
	testRow <- x.auto[1, , drop = FALSE]
	testRow[, "weight"] <- testRow[, "weight"] / 2

	print(testRow)
	print(head(x.auto))
	print(pintPerfect(testRow, x.auto, y.auto))
	print(pintBootstrap(testRow, x.auto, y.auto))
}

# main()