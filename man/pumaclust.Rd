\name{pumaClust}
\alias{pumaClust}
\title{Propagate probe-level uncertainty in model-based clustering on gene expression data}
\description{
  This function clusters gene expression using a Gaussian mixture model including probe-level 
  measurement error. The inputs are gene expression levels and the probe-level standard
  deviation associated with expression measurement for each gene on each chip. The outputs 
  is the clustering results.
}
\usage{
pumaClust(e=NULL, se=NULL, efile=NULL, sefile=NULL, 
      subset=NULL, gsnorm=FALSE, clusters, 
      iter.max=100, nstart=10, eps=1.0e-6, del0=0.01)
}
\arguments{
  \item{e}{ either a valid \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object, or a data frame containing the expression level for each gene on each chip. }
  \item{se}{ data frame containing the standard deviation of gene expression levels. }
  \item{efile}{ character, the name of the file which contains gene expression measurements. }
  \item{sefile}{ character, the name of the file which contains the standard deviation of gene expression measurements. }
  \item{subset}{ vector specifying the row number of genes which are clustered on.}
  \item{gsnorm}{ logical specifying whether do global scaling normalisation or not. }
  \item{clusters}{integer, the number of clusters. }
  \item{iter.max}{ integer, the maximum number of iterations allowed in the parameter initialisation. }
  \item{nstart}{ integer, the number of random sets chosen in the parameter initialisation. }
  \item{eps}{ numeric, optimisation parameter. }
  \item{del0}{ numeric, optimisation parameter. }
}
\details{
  The input data is specified either as an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object (in which case se, efile and sefile will be ignored), or by e and se, or by efile and sefile.
 
}
\value{
  The result is a list with components
  
  cluster: vector, containing the membership of clusters for each gene;
  centers: matrix, the center of each cluster;
  centersigs: matrix, the center variance of each cluster;
  likelipergene: matrix, the likelihood of belonging to each cluster for each gene;
  bic: numeric, the Bayesian Information Criterion score.
  
}
\references{
Liu,X., Lin,K.K., Andersen,B., and Rattray,M. (2006) 
Propagating probe-level uncertainty in model-based gene expression clustering, 
BMC Bioinformatics, 8(98).

Liu,X., Milo,M., Lawrence,N.D. and Rattray,M. (2005) A tractable probabilistic model for 
Affymetrix probe-level analysis across multiple chips, Bioinformatics, 21(18):3637-3644.
}
\author{ Xuejun Liu, Magnus Rattray }
\seealso{ Related method \code{\link{mmgmos}} and \code{\link[puma]{pumaClustii}}}
\examples{
  data(Clust.exampleE)
  data(Clust.exampleStd)
  pumaClust.example<-pumaClust(Clust.exampleE,Clust.exampleStd,clusters=7)
}
\keyword{ manip }
\keyword{ models }
