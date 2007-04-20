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
setMethod("plot",signature(x="pumaPCARes",y="missing"),
          function(x,y,...) {
              .plot.pumaPCARes(x,...)
          })

setClass("exprReslt"
,	representation(
		prcfive="matrix",
		prctwfive="matrix",
		prcfifty="matrix",
		prcsevfive="matrix",
		prcninfive="matrix"
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
if( !isGeneric("se.exprs") )
	setGeneric("se.exprs", function(object) standardGeneric("se.exprs"))
setMethod("se.exprs", "exprReslt", function(object) object@se.exprs)

if( !isGeneric("se.exprs<-") )
	setGeneric("se.exprs<-", function(object, value)
		standardGeneric("se.exprs<-"))

setReplaceMethod("se.exprs", "exprReslt",
                 function(object, value) {
                   object@se.exprs <- value
                   return(object)
                 })

##define a generic for obtaining the data
if( !isGeneric("prcfive") )
	setGeneric("prcfive", function(object) standardGeneric("prcfive"))
setMethod("prcfive", "exprReslt", function(object) object@prcfive)

if( !isGeneric("prcfive<-") )
	setGeneric("prcfive<-", function(object, value)
		standardGeneric("prcfive<-"))

setReplaceMethod("prcfive", "exprReslt",
                 function(object, value) {
                   object@prcfive <- value
                   return(object)
                 })

##define a generic for obtaining the data
if( !isGeneric("prctwfive") )
	setGeneric("prctwfive", function(object) standardGeneric("prctwfive"))
setMethod("prctwfive", "exprReslt", function(object) object@prctwfive)

if( !isGeneric("prctwfive<-") )
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
if( !isGeneric("prcfifty") )
	setGeneric("prcfifty", function(object) standardGeneric("prcfifty"))
setMethod("prcfifty", "exprReslt", function(object) object@prcfifty)

if( !isGeneric("prcfifty<-") )
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
if( !isGeneric("prcsevfive") )
	setGeneric("prcsevfive", function(object) standardGeneric("prcsevfive"))
setMethod("prcsevfive", "exprReslt", function(object) object@prcsevfive)

if( !isGeneric("prcsevfive<-") )
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
if( !isGeneric("prcninfive") )
	setGeneric("prcninfive", function(object) standardGeneric("prcninfive"))
setMethod("prcninfive", "exprReslt", function(object) object@prcninfive)

if( !isGeneric("prcninfive<-") )
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
	 
if( !isGeneric("write.reslts") )
	setGeneric("write.reslts", function(x,...) standardGeneric("write.reslts"))
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


