create_eset_r <- function(eset, r)
{
	numOfFactors <- numOfFactorsToUse(removeUninformativeFactors(eset))
	pData_eset <- as.data.frame(unique(pData(eset)[,1:numOfFactors]))
	rownames(pData_eset) <- colnames(r)[1:(dim(r)[2]/2)]
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
