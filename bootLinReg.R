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

pressReg <- function(xmat0, yvec, print.ls = TRUE, resid.plot = TRUE, quadChoice = FALSE) {
	# setup plotting area
	par(mfrow = c(1, 2))

	# run leaps
	library(leaps)
	xmat <- matrix.2ndorder.make(xmat0, only.quad = quadChoice)
	leaps.str <- leaps(xmat, yvec)

	# rank data by Cp (calculate Cp - p and order)
	z1 <- leaps.str$Cp - leaps.str$size
	o1 <- order(z1)
	parvec <- (leaps.str$which[o1, ])[1, ]
	xsubset <- xmat[, parvec]
	Cpval <- z1[o1][[1]]

	ls.str0 <- regpluspress(xsubset, yvec)
	if (print.ls) {
		ls.print(ls.str0)
	}

	npar <- sum(parvec)
	mpse <- ls.str0$press / (length(yvec) - (npar + 1))
	coef <- ls.str0$coef

	print(sprintf('Press = %f', ls.str0$press))
	print(sprintf('MPSE = %f', mpse))
	print(sprintf('Cp - p = %f', Cpval))

	MPSE <- floor(1000 * mpse / 1000)
	ypred <- xsubset %*% c(coef[-1]) + coef[1]

	plot(ypred, yvec, 
		main = sprintf("MPSE = %f, Cp - p = %f", 
						MPSE, Cpval)
	)

	if (resid.plot) {
		plot(ypred, ls.str0$resid)
	}
}