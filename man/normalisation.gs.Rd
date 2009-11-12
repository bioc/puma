\name{normalisation.gs}
\alias{normalisation.gs}
\title{ Global scaling normalisation }
\description{
This function is only included for backwards compatibility with the \pkg{pplr} package. This function is now superceded by \code{\link{pumaNormalize}}.

This function does the global scaling normalisation.
}
\usage{
normalisation.gs(x)
}
\arguments{
  \item{x}{ a matrix or data frame which contains gene expression level on log2 scale. }
}
\details{
  Each row of x is related to a gene and each column is related to a chip.
}
\value{
  The return matrix is in the same format as the input x.
}
\author{ Xuejun Liu, Marta Milo, Neil D. Lawrence, Magnus Rattray }
\seealso{ See Also as \code{\link{bcomb}} and \code{\link{hcomb}} }
\examples{
data(exampleE)
exampleE.normalised<-normalisation.gs(exampleE)
data(Clust.exampleE)
Clust.exampleE.normalised<-normalisation.gs(Clust.exampleE)
}
\keyword{ manip }
