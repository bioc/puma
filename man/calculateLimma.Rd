\name{calculateLimma}
\alias{calculateLimma}
\title{Calculate differential expression between conditions using limma}
\description{
Runs a default analysis using the \pkg{limma} package. Automatically creates design and contrast matrices if not specified. This function is useful for comparing \pkg{limma} results with those of \code{\link{pumaDE}}.
}
\usage{
calculateLimma(
	eset
,	design.matrix = createDesignMatrix(eset)
,	contrast.matrix = createContrastMatrix(eset)
,	direction="either"
)
}
\arguments{
	\item{eset}{ An object of class \code{\link[Biobase]{ExpressionSet}} }
	\item{design.matrix}{ A design matrix }
	\item{contrast.matrix}{ A contrast matrix }
	\item{direction}{ A string ("either", "up" or "down") }
}
\details{
The \code{eset} argument must be supplied, and must be a valid \code{\link[Biobase]{ExpressionSet}} object. Design and contrast matrices can be supplied, but if not, default matrices will be used. These should usually be sufficient for most analyses.
}
\value{
  A list with two elements. The first element (\code{p}) is a matrix of p-values. Each column represent one contrast. Within each column the p-values are ordered. The second element (\code{genes}) is a matrix of row numbers, which can be used to map p-values back to probe sets.
}
\author{ Richard D. Pearson }
\seealso{Related methods \code{\link{pumaDE}}, \code{\link{createDesignMatrix}} and \code{\link{createContrastMatrix}}}
\examples{
	library(puma)
	data(affybatch.example)
	pData(affybatch.example) <- data.frame("level"=c("twenty","twenty","ten")
	    , "batch"=c("A","B","A"), row.names=rownames(pData(affybatch.example)))
	eset_rma <- rma(affybatch.example)
	limmaRes <- calculateLimma(eset_rma)
	limmaRes$p[1:6]
	plotErrorBars(eset_rma, limmaRes$genes[1])
}
\keyword{manip}