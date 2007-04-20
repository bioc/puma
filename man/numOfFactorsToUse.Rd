\name{numOfFactorsToUse}
\alias{numOfFactorsToUse}
\title{Determine number of factors to use from an ExpressionSet}
\description{
This is really an internal function used to determine how many factors to use in design and contrast matrices}
\usage{
numOfFactorsToUse(eset)
}
\arguments{
	\item{eset}{ An object of class \code{\link[Biobase]{ExpressionSet}}. }
}
\value{
  An integer denoting the number of factors to be used.
}
\author{ Richard D. Pearson }
\seealso{Related methods \code{\link{createDesignMatrix}} and  \code{\link{createContrastMatrix}}}
\examples{
	data(affybatch.example)
	pData(affybatch.example) <- data.frame("level"=c("twenty","twenty","ten")
	    , "batch"=c("A","B","A"), row.names=rownames(pData(affybatch.example)))
	eset_mmgmos <- mmgmos(affybatch.example)
	numOfFactorsToUse(eset_mmgmos)
}
\keyword{manip}
