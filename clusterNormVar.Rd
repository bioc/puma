\name{clusterNormVar}
\alias{clusterNormVar}
\title{ Adjusting expression variance for zero-centered normalisation }
\description{
  This function adjusts the variance of the gene expression according to the zero-centered normalisation.
}
\usage{
clusterNormVar(x)
}
\arguments{
  \item{x}{ a vector which contains the variance of gene expression level on log2 scale. }
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
data(Clust.exampleStd)
Clust.exampleVar<-Clust.exampleStd^2
Clust.exampleStd.centered<-t(apply(cbind(Clust.exampleE,Clust.exampleVar), 1, clusterNormVar))
}
\keyword{ manip }
