# smoothers.R
# Runs a running mean, running median, and combines them on a given data set

testFunc <- function(x) {
	return( x^2 )
}

errorFunc <- function(x) {
	return( rnorm(1, mean=0, sd=500) )
}

gaussianNoiseOnFunc <- function(xlist, testFunc) {
	rawData <- lapply(xlist, testFunc)
	errorData <- lapply(xlist, errorFunc)
	finalData <- mapply(function(x, y) {
		return( x + y )
	}, rawData, errorData)
	return( finalData )
}

x0 <- c(1:100)
y0 <- gaussianNoiseOnFunc(x0, testFunc)

runningMean <- function(xlist, ylist, windowSize=10) {
	ylist <- c(rep(ylist[1], times = floor(windowSize / 2)), ylist, rep(ylist[length(ylist)], times = ceiling(windowSize / 2)))
	startInd <- floor(windowSize / 2)
	endInd <- length(ylist) + ceiling(windowSize / 2)
	ymean <- lapply(c(startInd:endInd))
	return( ymean )
}

runningMedian <- function(xlist, ylist) {
	return( ylist )
}

mySmooth <- function(xlist, ylist, windowSize=10, p=1) {
	ymean <- runningMean(xlist, ylist, windowSize)
	ymedian <- runningMedian(xlist, ylist, windowSize)

	ycomb <- p * ymean + (1 - p) * ymedian
	#
	#set up plotting region
	#
	par(mfrow=c(2,2))
	#
	#plot x0 vs y
	#
	plot(xlist,ylist)
	lines(x0,ycomb,col=2)
}