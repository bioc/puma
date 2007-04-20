calculateLimma <- function (
	eset
,	design.matrix = createDesignMatrix(eset)
,	contrast.matrix = createContrastMatrix(eset)
,	direction="either"
)
{
	numOfGenes <- length(featureNames(eset))
	numOfContrasts <- dim(contrast.matrix)[2]
	fit <- lmFit(eset, design.matrix)
	fit2 <- contrasts.fit(fit, contrast.matrix)
	fit2 <- eBayes(fit2)
	tt <- list()
	p <- matrix(0,numOfGenes, numOfContrasts)
	genes <- matrix(0,numOfGenes, numOfContrasts)
	for(i in 1:numOfContrasts)
	{
		tt[[i]] <- topTable(fit2, coef=i, number=numOfGenes)
		if(direction=="either")
			p[,i] <- tt[[i]]$P.Value
		else if(direction=="up")
			p[,i] <- ( tt[[i]]$P.Value * as.integer(tt[[i]]$M > 0) ) +  ( (1-tt[[i]]$P.Value) * as.integer(tt[[i]]$M < 0) )
		else if(direction=="down")
			p[,i] <- ( tt[[i]]$P.Value * as.integer(tt[[i]]$M < 0) ) +  ( (1-tt[[i]]$P.Value) * as.integer(tt[[i]]$M > 0) )
		genes[,i] <- as.integer(rownames(tt[[i]]))
	}
	list(p=p, genes=genes)
}
