calculateLimma <- function (
	eset
,	design.matrix = createDesignMatrix(eset)
,	contrast.matrix = createContrastMatrix(eset)
)
{
	require(limma)
	
	numOfGenes <- length(featureNames(eset))
	numOfContrasts <- dim(contrast.matrix)[2]
	fit <- lmFit(eset, design.matrix)
	fit2 <- contrasts.fit(fit, contrast.matrix)
	fit2 <- eBayes(fit2)
	p.value <- fit2$p.value
	fold.change <- fit2$coefficients
	new("DEResult", statistic=p.value, FC=fold.change
	, statisticDescription="limma p.value", DEMethod="calculateLimma")
}

calculateFC <- function (
	eset
,	design.matrix = createDesignMatrix(eset)
,	contrast.matrix = createContrastMatrix(eset)
)
{
	p <- matrix(0, dim(exprs(eset))[1], dim(contrast.matrix)[2]
	, dimnames=list(rownames(exprs(eset)), colnames(contrast.matrix)))
    fold.change <- matrix(0, dim(exprs(eset))[1], dim(contrast.matrix)[2]
	, dimnames=list(rownames(exprs(eset)), colnames(contrast.matrix)))
    for(i in 1:dim(contrast.matrix)[2])
	{
		plus_contrasts <- which(contrast.matrix[,i]==1)
		minus_contrasts <- which(contrast.matrix[,i]==-1)
		plus_arrays <- as.matrix(
			which(design.matrix[,plus_contrasts]==1,arr.ind=TRUE))[,1]
		minus_arrays <- as.matrix(
			which(design.matrix[,minus_contrasts]==1,arr.ind=TRUE))[,1]
		fold.change[,i] <-	rowMeans(as.matrix(exprs(eset)[,plus_arrays])) -
			rowMeans(as.matrix(exprs(eset)[,minus_arrays]))
	}
	p <- abs(fold.change)
	new("DEResult", statistic=p, FC=fold.change
	, statisticDescription="abs(fold change)", DEMethod="calculateFC")
}

calculateTtest <- function (
	eset
,	design.matrix = createDesignMatrix(eset)
,	contrast.matrix = createContrastMatrix(eset)
)
{
    p <- matrix(0, dim(exprs(eset))[1], dim(contrast.matrix)[2]
	, dimnames=list(rownames(exprs(eset)), colnames(contrast.matrix)))
    fold.change <- matrix(0, dim(exprs(eset))[1], dim(contrast.matrix)[2]
	, dimnames=list(rownames(exprs(eset)), colnames(contrast.matrix)))
    for(i in 1:dim(contrast.matrix)[2])
	{
		plus_contrasts <- which(contrast.matrix[,i]==1)
		minus_contrasts <- which(contrast.matrix[,i]==-1)
		plus_arrays <- as.matrix(which(design.matrix[,plus_contrasts]==1,arr.ind=TRUE))[,1]
		minus_arrays <- as.matrix(which(design.matrix[,minus_contrasts]==1,arr.ind=TRUE))[,1]
		p[,i] <- apply( exprs(eset), 1
		,	function(x) t.test(x[plus_arrays],x[minus_arrays])$p.value)
		fold.change[,i] <-	rowMeans(as.matrix(exprs(eset)[,plus_arrays])) -
			rowMeans(as.matrix(exprs(eset)[,minus_arrays]))
	}
	new("DEResult", statistic=p, FC=fold.change
	, statisticDescription="T-test p.value", DEMethod="calculateTtest")
}

