examineDensities <- function(datamat=Caravan, resultIndex=86) {
	for (i in (c(1:length(datamat)))) {
		if (i != resultIndex) {
			plotname <- paste("plot", i, ".png", sep = "")
			print(plotname)
			print(dev.cur())
			png(plotname)
			print(densityplot(as.formula(paste("~", colnames(datamat)[i], "|", colnames(datamat)[resultIndex])), 
						data = as.data.frame(datamat)))
			print(dev.cur())
			dev.off()
			# readline()
		}
	}
}

regOnSubset <- function(data=Caravan, resultInd, inputRange, linkFam="probit") {
	formulaStr <- paste(colnames(data)[inputRange], collapse = " + ")
	formulaStr <- paste(resultInd, "~", formulaStr)
	# print(formulaStr)
	return( glm(as.formula(formulaStr), family = binomial(link = linkFam), data = data) )
}

indivRegs <- function(data) {
	pvals <- NULL
	for (i in c(1:85)) {
		pvals[i] <- coef(summary(regOnSubset(data, "Purchase", i)))[, 4][2]
	}
	names(pvals) <- colnames(data)[c(1:85)]

	return( pvals )
}

getTrainResults <- function(actual, glmModel) {
	trainResults <- predict(glmModel, type = "response") > 0.5
	predicted <- as.factor(ifelse(trainResults > 0.5, "Yes", "No"))
	table(actual, predicted)
}
