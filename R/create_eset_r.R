create_eset_r <- function(eset, r, design.matrix=createDesignMatrix(eset))
{
	numOfFactors <- numOfFactorsToUse(removeUninformativeFactors(eset))
	pData_eset <- as.data.frame(
		as.data.frame(unique(pData(removeUninformativeFactors(eset))[,1:numOfFactors]))[
			apply(unique(design.matrix),2,function(x) which(x==1)),]
		)
	colnames(pData_eset) <- colnames(pData(removeUninformativeFactors(eset)))[1:numOfFactors]
	rownames(pData_eset) <- colnames(design.matrix)
	exprs <- as.matrix(r[,1:(dim(r)[2]/2)])
	se.exprs <- as.matrix(r[,((dim(r)[2]/2)+1):(dim(r)[2])])
	colnames(exprs) <- NULL
	colnames(se.exprs) <- NULL
	eset_r <- new(
		'ExpressionSet'
	,	exprs=exprs
	,	se.exprs=se.exprs
	,	phenoData = new(
			"AnnotatedDataFrame"
		,	data=pData_eset
		,	varMetadata=data.frame(
				labelDescription=varLabels(phenoData(eset))[1:numOfFactors]
			)
				# ,	varMetadata=data.frame(
		# 		matrix(
		# 			varMetadata(eset)[1:numOfFactors,]
		# 		,	dimnames=list(
		# 				rownames(varMetadata(eset))[1:1]
		# 			,	"labelDescription"
		# 			)
		# 		)
		# 	)
		)
	)
  # pData(eset_r) <- pData_eset
	return(eset_r)
}

