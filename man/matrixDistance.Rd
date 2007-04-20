\name{matrixDistance}
\alias{matrixDistance}
\title{Calculate distance between two matrices}
\description{
This calculates the mean Euclidean distance between the rows of two matrices. It is used in the function \code{\link{pumaPCA}}
}
\usage{
matrixDistance(
    matrixA
,   matrixB
)
}
\arguments{
	\item{matrixA}{ the first matrix }
	\item{matrixB}{ the second matrix }
}
\value{
	A numeric giving the mean distance
}
\author{ Richard D. Pearson }
\seealso{Related class \code{\link{pumaPCA}}}
\examples{
	show(matrixDistance(matrix(1,2,2),matrix(2,2,2)))
}
\keyword{manip}
