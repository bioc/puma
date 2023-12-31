\name{calculateLimma}
\alias{calculateLimma}
\title{Calculate differential expression between conditions using limma}
\description{
Runs a default analysis using the \pkg{limma} package. Automatically creates design and contrast matrices if not specified. This function is useful for comparing \pkg{limma} results with those of other differential expression (DE) methods such as \code{\link{pumaDE}}.
}
\usage{
calculateLimma(
	eset
,	design.matrix = createDesignMatrix(eset)
,	contrast.matrix = createContrastMatrix(eset)
)
}
\arguments{
	\item{eset}{ An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} }
	\item{design.matrix}{ A design matrix }
	\item{contrast.matrix}{ A contrast matrix }
}
\details{
The \code{eset} argument must be supplied, and must be a valid \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object. Design and contrast matrices can be supplied, but if not, default matrices will be used. These should usually be sufficient for most analyses.
}
\value{
  An object of class \code{\link{DEResult}}.
}
\author{ Richard D. Pearson }
\seealso{Related methods \code{\link{pumaDE}}, \code{\link{calculateTtest}}, \code{\link{calculateFC}}, \code{\link{createDesignMatrix}} and \code{\link{createContrastMatrix}} and class \code{\link{DEResult}}}
\examples{
#if (require(affydata)) {
#	data(Dilution)
#	eset_rma <- affy:::rma(Dilution)
#	#	Next line used so eset_rma only has information about the liver factor
#	#	The scanner factor will thus be ignored, and the two arrays of each level
#	#	of the liver factor will be treated as replicates
#	pData(eset_rma) <- pData(eset_rma)[,1, drop=FALSE]
#	limmaRes <- calculateLimma(eset_rma)
#	topGeneIDs(limmaRes,numberOfGenes=6)
#	plotErrorBars(eset_rma, topGenes(limmaRes))
#}
}
\keyword{manip}
