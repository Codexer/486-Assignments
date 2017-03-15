my.hat <- function(x) {
	x1 <- cbind(1, x)
	x1 %*% solve(t(x1) %*% x1) %*% t(x1)
}

regpluspress <- function(x, y) {
	ls.str <- lsfit(x, y)
	ls.str$leverage <- hat(x)
	ls.str$press <- sum( (ls.str$resid / (1 - ls.str$leverage))^2 )
	return( ls.str )
}

pressRegAdj <- function(xmat, yvec, ncheck = 10, print.ls = FALSE, resid.plot = TRUE) {
	# setup plotting area
	# n1 <- ceiling(sqrt(ncheck))
	par(mfrow = c(1, 2))

	# add hp^2, displacement^2
	xmat <- cbind(xmat, xmat[, "horsepower"]^2)
	xmat <- cbind(xmat, xmat[, "displacement"]^2)

	# remove origin
	xmat <- xmat[, -7]

	# run leaps
	library(leaps)
	leaps.str <- leaps(xmat, yvec)

	# rank data by Cp (calculate Cp - p and order)
	z1 <- leaps.str$Cp - leaps.str$size
	o1 <- order(z1)
	matwhich <- (leaps.str$which[o1, ])[1:ncheck, ]
	z2 <- z1[o1][1:ncheck]

	for (i in 1:ncheck) {
		ls.str0 <- regpluspress(xmat[, matwhich[i, ]], yvec)
		if (print.ls) {
			ls.print(ls.str0)
		}

		parvec <- matwhich[i, ]
		npar <- sum(parvec)
		mpse <- ls.str0$press / (length(yvec) - (npar + 1))

		print(sprintf('Trial %d', i))
		print(sprintf('Press = %f', ls.str0$press))
		print(sprintf('MPSE = %f', mpse))
		print(sprintf('Cp - p = %f', z2[i]))

		MPSE <- floor(1000 * mpse / 1000)
		xhat <- my.hat(xmat[, matwhich[i, ]])
		ypred <- xhat %*% yvec

		plot(ypred, yvec, 
			main = sprintf("i = %d, MPSE = %f, Cp - p = %f", 
							i, MPSE, z2[i])
		)
		if (resid.plot) {
			plot(ypred, ls.str0$resid)
		}
		
		readline()
	}
}

pressRegComp <- function(xmat, yvec, ncheck = 10, print.ls = TRUE, resid.plot = TRUE) {
	# setup plotting area
	# n1 <- ceiling(sqrt(ncheck))
	par(mfrow = c(1, 2))

	# run leaps
	library(leaps)
	leaps.str <- leaps(xmat, yvec)

	# rank data by Cp (calculate Cp - p and order)
	z1 <- leaps.str$Cp - leaps.str$size
	o1 <- order(z1)
	matwhich <- (leaps.str$which[o1, ])[1:ncheck, ]
	z2 <- z1[o1][1:ncheck]

	for (i in 1:ncheck) {
		parvec <- matwhich[i, ]
		whichX <- xmat[, parvec]
		numVars <- dim(whichX)[[2]]

		for (j in 1:numVars) {
			currX <- whichX[, j]		
			ls.str0 <- regpluspress(currX, yvec)
			if (print.ls) {
				ls.print(ls.str0)
			}

			npar <- sum(parvec)
			mpse <- ls.str0$press / (length(yvec) - (npar + 1))

			print(sprintf('Trial %d', i))
			print(sprintf('Press = %f', ls.str0$press))
			print(sprintf('MPSE = %f', mpse))
			print(sprintf('Cp - p = %f', z2[i]))
			print(colnames(whichX))

			MPSE <- floor(1000 * mpse / 1000)
			xhat <- my.hat(currX)
			ypred <- xhat %*% yvec

			plot(ypred, yvec, 
				main = sprintf("only %s", colnames(whichX)[[j]])
			)
			if (resid.plot) {
				plot(ypred, ls.str0$resid)
			}
			
			readline()
		}
	}
}