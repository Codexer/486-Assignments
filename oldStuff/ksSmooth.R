my.smooth.forKS <- function(x.str, xindex, yindex, ind.sqrt = TRUE, debug = FALSE) {
	if (debug) {
		par(mfrow = c(2, 2))
		plot(x.str[[xindex]], sqrt(x.str[[yindex]]), main = "Raw data and smooths")
		lines(smooth.spline(x.str[[xindex]], sqrt(x.str[[yindex]])))
	}

	if (ind.sqrt) {
		smsp.strcv <- smooth.spline(x.str[[xindex]], sqrt(x.str[[yindex]]))
		smspcv.resid <- sqrt(x.str[[yindex]]) - approx(smsp.strcv$x, smsp.strcv$y, x.str[[xindex]])$y
	} else {
		smsp.strcv <- smooth.spline(x.str[[xindex]], x.str[[yindex]])
		smspcv.resid <- x.str[[yindex]] - approx(smsp.strcv$x, smsp.strcv$y, x.str[[xindex]])$y
	}

	sd.resid <- sqrt(sum(smspcv.resid^2) / ( length(x.str[[1]]) - smsp.strcv$df ))
	stud.resid <- smspcv.resid / sd.resid
	ks.str <- ks.test(stud.resid, pnorm)$statistic
	D <- ks.str$statistic
	Pval <- ks.str$p.value
	my.smooth <- approx(smsp.strcv$x, smsp.strcv$y, x.str[[xindex]])$y
	list(D = D, raw.resid = smspcv.resid, sd.resid = sd.resid, smooth = my.smooth, P = P.val)
}

my.KS.pboot.normal <- function(mydata, x.index, y.index, nboot = 1000, mydf = NULL) {
	str0 <- my.smooth.forKS(mydata, x.index, y.index)
	Ddist <- NULL
	base.smooth <- str0$smooth
	base.sd <- str0$sd.resid
	D0 <- str0$D
	my.bootdata <- mydata
	n1 <- length(base.smooth)

	for (i in 1:nboot) {
		bres <- rnorm(n1, 0, base.sd)
		boot.dat <- ((base.smooth + bres))
		# print(boot.dat)
		my.bootdata[[y.index]] <- boot.dat
		bstr0 <- my.smooth.forKS(my.bootdata, x.index, y.index, FALSE)
		Ddist <- c(Ddist, bstr0$D)
	}

	Pval <- sum(Ddist > D0) / nboot
	Pval
}