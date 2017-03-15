regressOnNOAA <- function(tempcut=0.4) {
	NOAA <- read.csv('NOAA+GISS.csv')
	NOAA <- NOAA[order(NOAA[["delta.temp"]]), ]
	tempData <- NOAA[["delta.temp"]]
	disasterData <- sqrt(NOAA[["X.disaster"]])

	cutoff <- length(tempData[tempData <= tempcut])

	sec1 <- c(1:cutoff)
	sec2 <- c((cutoff + 1):length(tempData))

	nf <- layout(matrix(c(1,1,1,1,2,3), 3, 2, byrow = TRUE))

	plot(tempData, disasterData, xlab="Temperature", ylab="sqrt(Number of disasters)", main="Billion dollar disasters per year")

	line1 <- lsfit(tempData[sec1], disasterData[sec1], intercept=TRUE)
	coef1 <- line1$coef
	x0 <- tempData[[1]]
	y0 <- coef1[1] + (x0 * coef1[-1])
	x1 <- tempData[[cutoff]]
	y1 <- coef1[1] + (x1 * coef1[-1])
	segments(x0, y0, x1, y1, col = "red")
	ls.print(line1)

	line2 <- lsfit(tempData[sec2], disasterData[sec2], intercept=TRUE)
	coef2 <- line2$coef
	x0 <- tempData[[cutoff + 1]]
	y0 <- coef2[1] + (x0 * coef2[-1])
	x1 <- tempData[[length(tempData)]]
	y1 <- coef2[1] + (x1 * coef2[-1])
	segments(x0, y0, x1, y1, col = "blue")
	ls.print(line2)

	plot(sec1, line1$res)
	plot(sec2, line2$res)

	readline()
}

lapply(0.1 * c(2:8), regressOnNOAA)
# regressOnNOAA()

dev.off()