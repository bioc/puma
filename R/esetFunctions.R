# se.exprs <- function (eset) assayDataElement(eset,"se.exprs")
# 
pumaNormalize <- function (
	eset
,	arrayScale = c("median", "none", "mean", "meanlog")
,	probesetScale = c("none", "mean", "median")
,	probesetNormalisation = NULL
,	replicates = list(1:dim(exprs(eset))[2])
)
{
	if(!(
		is(eset, "ExpressionSet") ||
		is(eset, "exprReslt")
	))
		stop("eset is not a valid ExpressionSet object!")
	E <- exprs(eset)
	SE <- assayDataElement(eset,"se.exprs")
	chipm <- rep(0,dim(exprs(eset))[2])
	if (arrayScale[1]=="mean")
	{
		expr <- as.data.frame(2^E)

		chipm <- apply(expr,2,mean)
		chipm <- chipm/chipm[1]
	}
	else if (arrayScale[1]=="median")
	{
		expr <- as.data.frame(2^E)

		chipm <- apply(expr,2,median)
		chipm <- chipm/chipm[1]
	}
	else if (arrayScale[1]=="meanlog")
	{
		chipm <- apply(E,2,mean)
		chipm <- chipm-chipm[1]
	}
	# expr <- as.matrix(log2(expr))
	if (
		arrayScale[1]=="mean" |
		arrayScale[1]=="median"
	)
	{
		for (i in 1:dim(exprs(eset))[2])
		{
			exprs(eset)[,i] <- exprs(eset)[,i]-log2(chipm[i])
			if(is(eset, "exprReslt"))
			{
				prcfive(eset)[,i] <- prcfive(eset)[,i]-log2(chipm[i])
				prctwfive(eset)[,i] <- prctwfive(eset)[,i]-log2(chipm[i])
				prcfifty(eset)[,i] <- prcfifty(eset)[,i]-log2(chipm[i])
				prcsevfive(eset)[,i] <- prcsevfive(eset)[,i]-log2(chipm[i])
				prcninfive(eset)[,i] <- prcninfive(eset)[,i]-log2(chipm[i])
			}
		}
	}
	if (
		arrayScale[1]=="meanlog" |
		arrayScale[1]=="medianlog"
	)
	{
		for (i in 1:dim(exprs(eset))[2])
		{
			exprs(eset)[,i] <- exprs(eset)[,i]-chipm[i]
			if(is(eset, "exprReslt"))
			{
				prcfive(eset)[,i] <- prcfive(eset)[,i]-chipm[i]
				prctwfive(eset)[,i] <- prctwfive(eset)[,i]-chipm[i]
				prcfifty(eset)[,i] <- prcfifty(eset)[,i]-chipm[i]
				prcsevfive(eset)[,i] <- prcsevfive(eset)[,i]-chipm[i]
				prcninfive(eset)[,i] <- prcninfive(eset)[,i]-chipm[i]
			}
		}
	}

	# if(arrayScale[1] == "mean")
	# 	E <- log2(sweep(2^E, 2, apply(2^E, 2, mean)))
	# if(arrayScale[1] == "median")
	# 	E <- log2(sweep(2^E, 2, apply(2^E, 2, median)))
	# if(arrayScale[1] == "meanlog")
	# 	E <- sweep(E, 2, apply(E, 2, mean))
	# if(arrayScale[1] == "medianlog")
	# 	E <- sweep(E, 2, apply(E, 2, median))
	E <- exprs(eset)
	# if(!is.null(probesetScale))
	# {
		# E_normd <- matrix(nrow=dim(E)[1], ncol=dim(E)[2])
		if(probesetScale[1]=="mean")
		{
			for(repl in replicates)
				E[,repl] <- sweep(E[,repl], 1, apply(E[,repl], 1, mean))
		}
		if(probesetScale[1]=="median")
		{
			for(repl in replicates)
				E[,repl] <- sweep(E[,repl], 1, apply(E[,repl], 1, median))
		}
		exprs(eset) <- E
	# }
	if(!is.null(probesetNormalisation))
	{
		assayDataElement(eset,"se.exprs") <- sqrt(
			t(
				apply(
					cbind(
						exprs(eset)
					,	assayDataElement(eset,"se.exprs")^2
					)
				,	1
				,	clusterNormVar
				)
			)
		)
		exprs(eset) <- t(apply(exprs(eset), 1, clusterNormE))
	}
	# exprs(eset) <- E
	return(eset)
}


# normalize se.exprs or exprs
# normalize array-level or probeset-level
# normalize value or value and se of value
# normalize to mean or median
