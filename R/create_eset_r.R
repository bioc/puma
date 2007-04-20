create_eset_r <- function(eset, r)
{
	numOfFactors <- numOfFactorsToUse(eset)
	pData_eset <- as.data.frame(unique(pData(eset)[,1:numOfFactors]))
	rownames(pData_eset) <- colnames(r)[1:(dim(r)[2]/2)]
	eset_r <- new(
		'ExpressionSet'
	,	exprs=as.matrix(r[,1:(dim(r)[2]/2)])
	,	se.exprs=as.matrix(r[,((dim(r)[2]/2)+1):(dim(r)[2])])
	,	phenoData = new(
			"AnnotatedDataFrame"
		,	data=pData_eset
		,	varMetadata=data.frame(
				matrix(
					varMetadata(eset)[1:numOfFactors,]
				,	dimnames=list(
						rownames(varMetadata(eset))[1:1]
					,	"labelDescription"
					)
				)
			)
		)
	)
  # pData(eset_r) <- pData_eset
	return(eset_r)
}
