\name{removeUninformativeFactors}
\alias{removeUninformativeFactors}
\title{Remove uninformative factors from the phenotype data of an ExpressionSet}
\description{
This is really an internal function used to remove uninformative factors from the phenotype data. Uninformative factors here are defined as those which have the same value for all arrays in the \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.}
\usage{
	removeUninformativeFactors(eset)
}
\arguments{
	\item{eset}{ An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. }
}
\value{
  An \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object with the same data as the input, except for a new \code{phenoData} slot.
}
\author{ Richard D. Pearson }
\seealso{Related methods \code{\link{createDesignMatrix}} and  \code{\link{createContrastMatrix}}}
\examples{
	eset_test <- new("ExpressionSet", exprs=matrix(rnorm(400,8,2),100,4))
	pData(eset_test) <- data.frame("informativeFactor"=c("A", "A", "B", "B"), "uninformativeFactor"=c("X","X","X","X"))
	eset_test2 <- removeUninformativeFactors(eset_test)
	pData(eset_test)
	pData(eset_test2)
}
\keyword{manip}
