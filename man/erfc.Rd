\name{erfc}
\alias{erfc}
\title{ The complementary error function }
\description{
  This function calculates the complementary error function of an input x.
}
\usage{
erfc(x)
}
\arguments{
  \item{x}{ a numeric, the input. }
}
\details{
  erfc is implemented using the function qnorm.
}
\value{
 The return is a numeric.
}
\author{ Xuejun Liu }
\seealso{ \code{\link{qnorm}}}
\examples{
  erfc(0.5)
}
\keyword{ math }