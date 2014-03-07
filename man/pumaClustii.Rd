\name{pumaClustii}
\alias{pumaClustii}
\title{Propagate probe-level uncertainty in robust t mixture clustering on replicated gene expression data}
\description{
  This function clusters gene expression by including uncertainties of gene expression measurements from probe-level analysis models and 
replicate information into a robust t mixture clustering model. The inputs are gene expression levels and the probe-level standard
  deviation associated with expression measurement for each gene on each chip. The outputs 
  is the clustering results.
}
\usage{
pumaClustii(e=NULL, se=NULL, efile=NULL, sefile=NULL, 
      subset=NULL, gsnorm=FALSE, mincls, maxcls, conds, reps, verbose=FALSE,
      eps=1.0e-5, del0=0.01)
}
\arguments{
  \item{e}{ data frame containing the expression level for each gene on each chip. }
  \item{se}{ data frame containing the standard deviation of gene expression levels. }
  \item{efile}{ character, the name of the file which contains gene expression measurements. }
  \item{sefile}{ character, the name of the file which contains the standard deviation of gene expression measurements. }
  \item{subset}{ vector specifying the row number of genes which are clustered on.}
  \item{gsnorm}{ logical specifying whether do global scaling normalisation or not. }
  \item{mincls}{integer, the minimum number of clusters. }
  \item{maxcls}{ integer, the maximum number of clusters. }
  \item{conds}{ integer, the number of conditions. }
  \item{reps}{ vector, specifying which condition each column of the input data matrix belongs to. }
  \item{verbose}{ logical value. If 'TRUE' messages about the progress of the function is printed. }
  \item{eps}{ numeric, optimisation parameter. }
  \item{del0}{ numeric, optimisation parameter. }
}
\details{
  The input data is specified either by e and se, or by efile and sefile.
 
}
\value{
  The result is a list with components
  
  cluster: vector, containing the membership of clusters for each gene;
  centers: matrix, the center of each cluster;
  centersigs: matrix, the center variance of each cluster;
  likelipergene: matrix, the likelihood of belonging to each cluster for each gene;
  optK: numeric, the optimal number of clusters.
  optF: numeric, the maximised value of target function.
  
}
\references{
Liu,X. and Rattray,M. (2009) 
Including probe-level measurement error in robust mixture clustering of replicated microarray gene expression, 
Statistical Application in Genetics and Molecular Biology, 9(1), Article 42.

Liu,X., Lin,K.K., Andersen,B., and Rattray,M. (2007) 
Propagating probe-level uncertainty in model-based gene expression clustering, 
BMC Bioinformatics, 8:98.

Liu,X., Milo,M., Lawrence,N.D. and Rattray,M. (2005) A tractable probabilistic model for 
Affymetrix probe-level analysis across multiple chips, Bioinformatics, 21(18):3637-3644.
}
\author{ Xuejun Liu }
\seealso{ Related method \code{\link[puma]{mmgmos}} and \code{\link[puma]{pumaclust}}}
\examples{
  data(Clustii.exampleE)
  data(Clustii.exampleStd)
  r<-vector(mode="integer",0)
  for (i in c(1:20))
    for (j in c(1:4))
      r<-c(r,i)
  cl<-pumaClustii(Clustii.exampleE,Clustii.exampleStd,mincls=6,maxcls=6,conds=20,reps=r,eps=1e-3)
}
\keyword{ manip }
\keyword{ models }
