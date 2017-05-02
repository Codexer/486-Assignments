interestPval <- function(pvals, Q = 0.05, debug.plot = T) {
	# needs some tests to work
	if (is.null(pvals) || is.na(pvals)) {
		print("No P values given. Aborting.")
		return( NULL )
	}

	# save x values
	m <- length(pvals)
	inds <- c(1:m)
	# calculate straight line based off Q parameter
	lineQ <- Q * (inds / m)
	# sort the pvals
	sortedVals <- sort(pvals)

	# find the largest below the line
	pCandidates <- (sortedVals < lineQ)
	largestInd <- (m + 1) - match(T, rev(pCandidates))

	if (debug.plot) {
		par(mfrow = c(1, 1))

		# draw every point as a faded black, slightly smaller, filled point
		plot(inds, sortedVals, main = "P values vs. Q line",
			col = rgb(0, 0, 0, 0.1), cex = 0.7, pch = 19)
		# draw the Q line with a blue color
		lines(inds, lineQ, 
			col = "blue")

		# special stuff gets drawn only if stuff actually exists
		if (!is.na(largestInd)) {
			# get all points below the line
			belowLine <- sortedVals[1:largestInd]
			xBelow <- c(1:largestInd)

			# draw every point below the line as a strong, red point
			points(xBelow, belowLine, 
				col = "red", cex = 0.7, pch = 19)
			# draw a vertical red line at the max val below the line
			abline(v = largestInd, 
				col = "red")
		}
	}

	# if there is nothing below the line, then abort
	if (is.na(largestInd)) {
		print("No P value found below the line. Aborting.")
		return( NULL )
	}

	# return the x value of every test that is significant
	return( inds[pvals <= sortedVals[largestInd]] )
}