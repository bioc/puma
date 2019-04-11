setClass("pumaPCAModel"
,	representation(
		sigma	=	"numeric"
    ,	m		=	"matrix"
    ,	Cinv	=	"matrix"
    ,	W		=	"matrix"
	,	mu		=	"matrix"
	)
,	prototype(
		sigma	=	0
    ,	m		=	matrix()
    ,	Cinv	=	matrix()
    ,	W		=	matrix()
	,	mu		=	matrix()
	)
)
# setAs("NULL", "pumaPCAModel", function(from, to) new(to))

setClass("pumaPCAExpectations"
,	representation(
		x			=	"matrix"
	,	xxT			=	"array"
	,	logDetCov	=	"numeric"
	)
,	prototype(
		x			=	matrix()
	,	xxT			=	array()
	,	logDetCov	=	0
	)
)
# setAs("NULL", "pumaPCAExpectations", function(from, to) new(to))

setClass("pumaPCARes"
,	representation(
	    model               = "pumaPCAModel"
    ,   expectations        = "pumaPCAExpectations"
    ,   varY                = "matrix"
    ,   Y                   = "matrix"
	,	phenoData			= "AnnotatedDataFrame"
    ,   timeToCompute       = "numeric"
    ,   numberOfIterations  = "numeric"
    ,   likelihoodHistory   = "list"
    ,   timingHistory       = "list"
    ,   modelHistory        = "list"
    ,   exitReason          = "character"
	)
,	prototype(
	    model               = new("pumaPCAModel")
    ,   expectations        = new("pumaPCAExpectations")
    ,   varY                = matrix()
    ,   Y                   = matrix()
	,	phenoData			= new("AnnotatedDataFrame")
    ,   timeToCompute       = 0
    ,   numberOfIterations  = 0
    ,   likelihoodHistory   = list()
    ,   timingHistory       = list()
    ,   modelHistory        = list()
    ,   exitReason          = ""
	)
)
# setAs("NULL", "pumaPCARes", function(from, to) new(to))
# setMethod("plot",signature(x="pumaPCARes",y="missing"),
#           function(x,y,...) {
#               .plot.pumaPCARes(x,...)
#           })

setGeneric("write.reslts", function(x,...) standardGeneric("write.reslts"))

setMethod(
	"write.reslts"
,	signature(x="pumaPCARes")
,	function(
		x
	,	file = "tmp"
	,	append = FALSE
	,	quote = FALSE
	,	sep = ","
	,	eol = "\n"
	,	na = "NA"
	,	dec = "."
	,	row.names = TRUE
	,	col.names = NA
	,	qmethod = c("escape", "double")
	)
	{
		write.table(
			x@model@W
		,	file = paste(file,".csv",sep="")
		,	append = append
		,	quote = quote
		,	sep = sep
		,	eol = eol
		,	na = na
		,	dec = dec
		,	row.names = rownames(pData(x@phenoData))
		,	col.names = col.names
		,	qmethod = qmethod
		)
	}
)

setClass("DEResult"
,	representation(
		statistic="matrix"
	,	FC="matrix"
	,	statisticDescription="character"
	,	DEMethod="character"
	)
,	prototype(
		statistic=matrix()
	,	FC=matrix()
	,	statisticDescription="unknown"
	,	DEMethod="unknown"
	)
)

setMethod("show", "DEResult",
          function(object) {
            cat("DEResult object:\n")
            cat("  DEMethod = ", object@DEMethod, "\n", sep="")
            cat("  statisticDescription = ", object@statisticDescription, "\n", sep="")
            cat("  statistic =", nrow(object@statistic), "probesets x", ncol(object@statistic), "contrasts\n")
          })

# if( !isGeneric("statistic"))
  setGeneric("statistic", function(object) standardGeneric("statistic"))

setMethod("statistic", "DEResult", function(object) object@statistic)

# if( !isGeneric("statistic<-"))
  setGeneric("statistic<-", function(object, value) standardGeneric("statistic<-"))

setReplaceMethod("statistic", signature=c("DEResult", "matrix"),
                 function(object, value) {
                   object@statistic <- as.matrix(value)
					return(object)
                 })

# if( !isGeneric("FC"))
  setGeneric("FC", function(object) standardGeneric("FC"))

setMethod("FC", "DEResult", function(object) object@FC)

# if( !isGeneric("FC<-"))
  setGeneric("FC<-", function(object, value) standardGeneric("FC<-"))

setReplaceMethod("FC", signature=c("DEResult", "matrix"),
                 function(object, value) {
                   object@FC <- as.matrix(value)
					return(object)
                 })


# if( !isGeneric("statisticDescription"))
  setGeneric("statisticDescription", function(object) standardGeneric("statisticDescription"))

setMethod("statisticDescription", "DEResult", function(object) object@statisticDescription)

# if( !isGeneric("statisticDescription<-"))
  setGeneric("statisticDescription<-", function(object, value) standardGeneric("statisticDescription<-"))

setReplaceMethod("statisticDescription", signature=c("DEResult", "character"),
                 function(object, value) {
                   object@statisticDescription <- as.character(value)
					return(object)
                 })

# if( !isGeneric("DEMethod"))
  setGeneric("DEMethod", function(object) standardGeneric("DEMethod"))

setMethod("DEMethod", "DEResult", function(object) object@DEMethod)

# if( !isGeneric("DEMethod<-"))
  setGeneric("DEMethod<-", function(object, value) standardGeneric("DEMethod<-"))

setReplaceMethod("DEMethod", signature=c("DEResult", "character"),
                 function(object, value) {
                   object@DEMethod <- as.character(value)
					return(object)
                 })



# if( !isGeneric("pLikeValues") )
	setGeneric("pLikeValues"
	,	function(object, contrast=1, direction="either") standardGeneric("pLikeValues"))
setMethod("pLikeValues", "DEResult",
	function(object, contrast=1, direction="either") {
		if(object@DEMethod=="pumaDE")
		{
			if(direction=="either")
				return(1-abs(2*(object@statistic[,contrast]-0.5)))
			if(direction=="up")
				return(1-object@statistic[,contrast])
			if(direction=="down")
				return(object@statistic[,contrast])
		}
		if(object@DEMethod=="calculateLimma" || object@DEMethod=="calculateTtest")
		{
			if(direction=="either")
				return(object@statistic[,contrast])
			if(direction=="up")
				return( ( ( (object@statistic[,contrast]-1) * sign(object@FC[,contrast]) )
					/ 2 ) + 0.5)
			if(direction=="down")
				return( ( ( (1-object@statistic[,contrast]) * sign(object@FC[,contrast]) )
				 	/ 2 ) + 0.5)
		}
		if(object@DEMethod=="calculateFC")
		{
			if(direction=="either")
				return(1-(rank(object@statistic[,contrast])-1)
					/	(length(object@statistic[,contrast])-1))
			if(direction=="up")
				return(1-(rank(object@FC[,contrast])-1)
					/	(length(object@FC[,contrast])-1))
			if(direction=="down")
				return((rank(object@FC[,contrast])-1)
					/	(length(object@FC[,contrast])-1))
		}
		if(object@DEMethod=="calculateCyberT")
		{
			if(direction=="either")
				return(1-(rank(abs(object@statistic[,contrast]))-1)
					/	(length(object@statistic[,contrast])-1))
			if(direction=="up")
				return((rank(object@statistic[,contrast])-1)
					/	(length(object@statistic[,contrast])-1))
			if(direction=="down")
				return(1-(rank(object@statistic[,contrast])-1)
					/	(length(object@statistic[,contrast])-1))
		}
     })

# if( !isGeneric("topGenes") )
	setGeneric("topGenes", function(object, numberOfGenes=1, contrast=1, direction="either") standardGeneric("topGenes"))
setMethod("topGenes", "DEResult",
	function(object, numberOfGenes=1, contrast=1, direction="either") {
		sort(pLikeValues(object, contrast, direction), method="quick"
		, index.return=TRUE)$ix[1:numberOfGenes]
	})

# if( !isGeneric("topGeneIDs") )
	setGeneric("topGeneIDs", function(object, numberOfGenes=1, contrast=1, direction="either") standardGeneric("topGeneIDs"))
setMethod("topGeneIDs", "DEResult",
	function(object, numberOfGenes=1, contrast=1, direction="either") {
		rownames(object@statistic)[sort(pLikeValues(object, contrast, direction), method="quick"
		, index.return=TRUE)$ix[1:numberOfGenes]]
	})

# if( !isGeneric("numberOfProbesets") )
	setGeneric("numberOfProbesets", function(object) standardGeneric("numberOfProbesets"))
setMethod("numberOfProbesets", "DEResult",
	function(object) {
		nrow(object@statistic)
	})

# if( !isGeneric("numberOfGenes") )
	setGeneric("numberOfGenes", function(object) standardGeneric("numberOfGenes"))
setMethod("numberOfGenes", "DEResult",
	function(object) {
		nrow(object@statistic)
	})

# if( !isGeneric("numberOfContrasts") )
	setGeneric("numberOfContrasts", function(object) standardGeneric("numberOfContrasts"))
setMethod("numberOfContrasts", "DEResult",
	function(object) {
		ncol(object@statistic)
	})

# if( !isGeneric("write.reslts") )
	# setGeneric("write.reslts", function(x,...) standardGeneric("write.reslts"))
setMethod(
	"write.reslts"
,	signature(x="DEResult")
,	function(
		x
	,	file = "tmp"
	,	append = FALSE
	,	quote = FALSE
	,	sep = ","
	,	eol = "\n"
	,	na = "NA"
	,	dec = "."
	,	row.names = TRUE
	,	col.names = NA
	,	qmethod = c("escape", "double")
	)
	{
		write.table(
			statistic(x)
		,	file = paste(file,"_statistics.csv",sep="")
		,	append = append
		,	quote = quote
		,	sep = sep
		,	eol = eol
		,	na = na
		,	dec = dec
		,	row.names = row.names
		,	col.names = col.names
		,	qmethod = qmethod
		)
		write.table(
			FC(x)
		,	file = paste(file,"_FCs.csv",sep="")
		,	append = append
		,	quote = quote
		,	sep = sep
		,	eol = eol
		,	na = na
		,	dec = dec
		,	row.names = row.names
		,	col.names = col.names
		,	qmethod = qmethod
		)
	}
)

setClass("exprReslt"
,	representation(
		prcfive="matrix",
		prctwfive="matrix",
		prcfifty="matrix",
		prcsevfive="matrix",
		prcninfive="matrix",
		exprs="matrix",   #this
		se.exprs="matrix",
		description="MIAME",
		annotation="character",
		notes="character",
		cdfName="character",
		nrow="numeric",
		ncol="numeric"
	)
,	prototype=list(
		prcfive=matrix(nr=0,nc=0),
		prctwfive=matrix(nr=0,nc=0),
		prcfifty=matrix(nr=0,nc=0),
		prcsevfive=matrix(nr=0,nc=0),
		prcninfive=matrix(nr=0,nc=0),
		exprs=matrix(nr=0,nc=0),
		se.exprs = matrix(nr=0,nc=0),
		description=new("MIAME"),
		annotation="",
		notes="",
		cdfName="",
		nrow=0,
		ncol=0
	)
,	contains="ExpressionSet"
)
	 
##define a generic for obtaining the data
# if( !isGeneric("se.exprs") )
#	setGeneric("se.exprs", function(object) standardGeneric("se.exprs"))
setMethod("se.exprs", "exprReslt", function(object) assayDataElement(object,"se.exprs"))

# if( !isGeneric("se.exprs<-") )
#	setGeneric("se.exprs<-", function(object, value)
#		standardGeneric("se.exprs<-"))

setReplaceMethod("se.exprs", "exprReslt",
                 function(object, value) {
                   assayDataElement(object,"se.exprs") <- value
                   return(object)
                 })

##define a generic for obtaining the data
# if( !isGeneric("prcfive") )
	setGeneric("prcfive", function(object) standardGeneric("prcfive"))
setMethod("prcfive", "exprReslt", function(object) object@prcfive)

# if( !isGeneric("prcfive<-") )
	setGeneric("prcfive<-", function(object, value)
		standardGeneric("prcfive<-"))

setReplaceMethod("prcfive", "exprReslt",
                 function(object, value) {
                   object@prcfive <- value
                   return(object)
                 })

##define a generic for obtaining the data
# if( !isGeneric("prctwfive") )
	setGeneric("prctwfive", function(object) standardGeneric("prctwfive"))
setMethod("prctwfive", "exprReslt", function(object) object@prctwfive)

# if( !isGeneric("prctwfive<-") )
	setGeneric("prctwfive<-", function(object, value)
    	standardGeneric("prctwfive<-"))

setReplaceMethod(
	"prctwfive"
,	"exprReslt"
,	function(object, value) {
		object@prctwfive <- value
		return(object)
	}
)
  ##define a generic for obtaining the data
# if( !isGeneric("prcfifty") )
	setGeneric("prcfifty", function(object) standardGeneric("prcfifty"))
setMethod("prcfifty", "exprReslt", function(object) object@prcfifty)

# if( !isGeneric("prcfifty<-") )
	setGeneric("prcfifty<-", function(object, value)
		standardGeneric("prcfifty<-"))

setReplaceMethod(
	"prcfifty"
,	"exprReslt"
,	function(object, value) {
		object@prcfifty <- value
		return(object)
	}
)
  ##define a generic for obtaining the data
# if( !isGeneric("prcsevfive") )
	setGeneric("prcsevfive", function(object) standardGeneric("prcsevfive"))
setMethod("prcsevfive", "exprReslt", function(object) object@prcsevfive)

# if( !isGeneric("prcsevfive<-") )
	setGeneric("prcsevfive<-", function(object, value)
		standardGeneric("prcsevfive<-"))

setReplaceMethod(
	"prcsevfive"
,	"exprReslt"
,	function(object, value) {
		object@prcsevfive <- value
		return(object)
	}
)
  ##define a generic for obtaining the data
# if( !isGeneric("prcninfive") )
	setGeneric("prcninfive", function(object) standardGeneric("prcninfive"))
setMethod("prcninfive", "exprReslt", function(object) object@prcninfive)

# if( !isGeneric("prcninfive<-") )
	setGeneric("prcninfive<-", function(object, value)
		standardGeneric("prcninfive<-"))

setReplaceMethod(
	"prcninfive"
,	"exprReslt"
,	function(object, value) {
		object@prcninfive <- value
		return(object)
	}
)

setMethod("show", "exprReslt", function(object ) {
	dm <-dim(exprs(object))
	ngenes <- dm[1]
	nsamples <- dm[2]
	cat("Expression Set (exprReslt) with \n\t", ngenes, " genes\n\t", sep="")
	cat(nsamples, "samples\n\t")
	show(phenoData(object))
})
	 
# if( !isGeneric("write.reslts") )
	# setGeneric("write.reslts", function(x,...) standardGeneric("write.reslts"))
setMethod(
	"write.reslts"
,	signature(x="exprReslt")
,	function(
		x
	,	file = "tmp"
	,	append = FALSE
	,	quote = FALSE
	,	sep = ","
	,	eol = "\n"
	,	na = "NA"
	,	dec = "."
	,	row.names = TRUE
	,	col.names = NA
	,	qmethod = c("escape", "double")
	)
	{
		write.table(
			exprs(x)
		,	file = paste(file,"_exprs.csv",sep="")
		,	append = append
		,	quote = quote
		,	sep = sep
		,	eol = eol
		,	na = na
		,	dec = dec
		,	row.names = row.names
		,	col.names = col.names
		,	qmethod = qmethod
		)
		write.table(
			se.exprs(x)
		,	file = paste(file,"_se.csv",sep="")
		,	append = append
		,	quote = quote
		,	sep = sep
		,	eol = eol
		,	na = na
		,	dec = dec
		,	row.names = row.names
		,	col.names = col.names
		,	qmethod = qmethod
		)
		write.table(
			prcfive(x)
		,	file = paste(file,"_prctile5.csv",sep="")
		,	append = append
		,	quote = quote
		,	sep = sep
		,	eol = eol
		,	na = na
		,	dec = dec
		,	row.names = row.names
		,	col.names = col.names
		,	qmethod = qmethod
		)
		write.table(
			prctwfive(x)
		,	file = paste(file,"_prctile25.csv",sep="")
		,	append = append
		,	quote = quote
		,	sep = sep
		,	eol = eol
		,	na = na
		,	dec = dec
		,	row.names = row.names
		,	col.names = col.names
		,	qmethod = qmethod
		)
		write.table(
			prcfifty(x)
		,	file = paste(file,"_prctile50.csv",sep="")
		,	append = append
		,	quote = quote
		,	sep = sep
		,	eol = eol
		,	na = na
		,	dec = dec
		,	row.names = row.names
		,	col.names = col.names
		,	qmethod = qmethod
		)
		write.table(
			prcsevfive(x)
		,	file = paste(file,"_prctile75.csv",sep="")
		,	append = append
		,	quote = quote
		,	sep = sep
		,	eol = eol
		,	na = na
		,	dec = dec
		,	row.names = row.names
		,	col.names = col.names
		,	qmethod = qmethod
		)
		write.table(
			prcninfive(x)
		,	file = paste(file,"_prctile95.csv",sep="")
		,	append = append
		,	quote = quote
		,	sep = sep
		,	eol = eol
		,	na = na
		,	dec = dec
		,	row.names = row.names
		,	col.names = col.names
		,	qmethod = qmethod
		)
	}
)

setMethod(
	"write.reslts"
,	signature(x="ExpressionSet")
,	function(
		x
	,	file = "tmp"
	,	append = FALSE
	,	quote = FALSE
	,	sep = ","
	,	eol = "\n"
	,	na = "NA"
	,	dec = "."
	,	row.names = TRUE
	,	col.names = NA
	,	qmethod = c("escape", "double")
	)
	{
		write.table(
			exprs(x)
		,	file = paste(file,"_exprs.csv",sep="")
		,	append = append
		,	quote = quote
		,	sep = sep
		,	eol = eol
		,	na = na
		,	dec = dec
		,	row.names = row.names
		,	col.names = col.names
		,	qmethod = qmethod
		)
	}
)
