\name{calculateTtest}
\alias{calculateTtest}
\title{Calculate differential expression between conditions using T-test}
\description{
Automatically creates design and contrast matrices if not specified. This function is useful for comparing T-test results with those of other differential expression (DE) methods such as \code{\link{pumaDE}}.
}
\usage{
calculateTtest(
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
\seealso{Related methods \code{\link{pumaDE}}, \code{\link{calculateLimma}}, \code{\link{calculateFC}}, \code{\link{createDesignMatrix}} and \code{\link{createContrastMatrix}} and class \code{\link{DEResult}}}
\examples{
	eset_test <- new("ExpressionSet", exprs=matrix(rnorm(400,8,2),100,4))
	pData(eset_test) <- data.frame("class"=c("A", "A", "B", "B"))
	TtestRes <- calculateTtest(eset_test)
	plotErrorBars(eset_test, topGenes(TtestRes))
}
\keyword{manip}