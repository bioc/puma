\name{clusterNormE}
\alias{clusterNormE}
\title{ Zero-centered normalisation }
\description{
  This function normalise the data vector to have zero mean.
}
\usage{
clusterNormE(x)
}
\arguments{
  \item{x}{ a vector which contains gene expression level on log2 scale. }
}
\details{
  Vector x is related to a gene and each element is related to a chip.
}
\value{
  The return vector is in the same format as the input x.
}
\author{ Xuejun Liu,  Magnus Rattray }
\seealso{ See Also as \code{\link{pumaClust}} and \code{\link{pumaClustii}} }
\examples{
data(Clust.exampleE)
Clust.exampleE.centered<-t(apply(Clust.exampleE, 1, clusterNormE))
}
\keyword{ manip }
