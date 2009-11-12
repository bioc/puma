\name{numOfFactorsToUse}
\alias{numOfFactorsToUse}
\title{Determine number of factors to use from an ExpressionSet}
\description{
This is really an internal function used to determine how many factors to use in design and contrast matrices}
\usage{
numOfFactorsToUse(eset)
}
\arguments{
	\item{eset}{ An object of class
	  \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. } 
}
\value{
  An integer denoting the number of factors to be used.
}
\author{ Richard D. Pearson }
\seealso{Related methods \code{\link{createDesignMatrix}} and  \code{\link{createContrastMatrix}}}
\examples{
	#	Next 4 lines commented out to save time in package checks, and saved version used
    # if (require(affydata)) {
	#	data(Dilution)
	#	eset_mmgmos <- mmgmos(Dilution)
	# }
	data(eset_mmgmos)
	numOfFactorsToUse(eset_mmgmos)
}
\keyword{manip}
