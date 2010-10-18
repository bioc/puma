\name{pumaDEUnsorted}
\alias{pumaDEUnsorted}
\title{Return an unsorted matrix of PPLR values}
\description{
Returns the output from \code{\link{pumaDE}} as an unsorted matrix (i.e. sorted according to the original sorting in the ExpressionSet)
}
\usage{
pumaDEUnsorted(pp)
}
\arguments{
	\item{pp}{ A list as output by \code{\link{pumaDE}}. }
}
\value{
  A matrix of PPLR values
}
\author{ Richard D. Pearson }
\seealso{Related method \code{\link{pumaDE}}}
\keyword{manip}
