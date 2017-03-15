# splines.R
# HW 1, doing various splines and F tests and fun stuff
# Matthew Tantoy (RUID: 163006720, NetID: mat322) for Stat 486

# Usage: splinesMain(frameName, indName1, indName2) where
# -- frameName is a data frame
# -- indName(1/2) are indices s.t. frameName[[indName(1/2)]] gets a column of the data frame

plotSpline <- function(dataFrame, index1, index2, colVal, dfVal = -1) {
	# grab the two columns of the data frame
	x <- dataFrame[[index1]]
	y <- dataFrame[[index2]]

	# use custom degrees of freedom if it's actually given
	if (dfVal == -1) {
		splineRet <- smooth.spline(x, y)
		titleStr <- "Smooth spline with default DF"
	} else {
		splineRet <- smooth.spline(x, y, df = dfVal)
		titleStr <- sprintf("Smooth spline with DF = %d", dfVal)
	}

	# plot the smooth fit on top of data points
	plot(x, y, xlab = index1, ylab = index2)
	title(main = titleStr)
	lines(splineRet, col = colVal)

	# calculate residuals using linear interpolation(?)
	splineRet$resid <- y - approx(splineRet$x, splineRet$y, x)$y

	return( splineRet )
}

# calculates F score using residuals
printFTest <- function(spline1, spline2, x) {
	SSF <- sum(spline1$resid^2)
	SSN <- sum(spline2$resid^2)

	Fscore <- ((SSN / SSF) / (spline1$df - 2)) / (SSF / length(x))

	# print F score and P value
	print(sprintf('F statistic: %f', Fscore))
	print(sprintf('P = %f', pf(Fscore, 1, length(x))))
}

splinesMain <- function(dataFrame, index1, index2) {
	par(mfrow = c(2, 2))

	spline1 <- plotSpline(dataFrame, index1, index2, colVal = "red")
	spline2 <- plotSpline(dataFrame, index1, index2, colVal = "blue", dfVal = 2)

	# quantile plots
	qqnorm(spline1$resid, main = "Normal Q-Q for default DF")
	qqnorm(spline2$resid, main = sprintf("Normal Q-Q for DF = %d", 2))

	printFTest(spline1, spline2, dataFrame[[index1]])
}