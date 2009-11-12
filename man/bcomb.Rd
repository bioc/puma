\name{bcomb}
\alias{bcomb}
\title{Combining replicates for each condition}
\description{
This function calculates the combined signal for each condition from replicates  using Bayesian models. The inputs are gene expression levels and the probe-level standard deviation associated with expression measurement for each gene on each chip. The outputs include gene expression levels and standard deviation for each condition.
This function was originally part of the \pkg{pplr} package. Although this function can be called directly, it is recommended to use the \code{\link[puma]{pumaComb}} function instead, which can work directly on \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} objects, and can automatically determine which arrays are replicates.
}
\usage{
bcomb(e, se, replicates, method=c("map","em"), 
      gsnorm=FALSE, nsample=1000, eps=1.0e-6)
}
\arguments{
  \item{e}{ a data frame containing the expression level for each gene on each chip. }
  \item{se}{ a data frame containing the standard deviation of gene expression levels. }
  \item{replicates}{ a vector indicating which chip belongs to which condition. }
  \item{method}{ character specifying the method algorithm used.}
  \item{gsnorm}{ logical specifying whether do global scaling normalisation or not. }
  \item{nsample}{ integer. The number of sampling in parameter estimation. }
  \item{eps}{ a numeric, optimisation parameter. }
}
\details{
  Each element in replicate represents the condition of the chip which is in the same column order as in the expression and standard deviation matrix files.
  
  Method "map" uses MAP of a hierarchical Bayesion model with Gamma prior on the between-replicate variance (Gelman et.al. p.285) and shares the same variance across conditions. This method is fast and suitable for the case where there are many conditions.
  
  Method "em" uses variational inference of the same hierarchical Bayesion model as in method "map" but with conjugate prior on between-replicate variance and shares the variance across conditions.
  
  The parameter nsample should be large enough to ensure stable parameter estimates. Should be at least 1000.
}
\value{
  The result is a data frame with components named 'M1', 'M2', and so on, which represent the mean expression values for condition 1, condition 2, and so on. It also has components named 'Std1', 'Std2', and so on, which represent the standard deviation of the gene expression values for condition 1, condtion 2, and so on.
}
\references{Gelman,A., Carlin,J.B., Stern,H.S., Rubin,D.B., Bayesian data analysis. London: Chapman & Hall; 1995.

Liu,X., Milo,M., Lawrence,N.D. and Rattray,M. (2006) Probe-level variances improve accuracy in detecting differential
gene expression, Bioinformatics, 22:2107-2113.}
\author{ Xuejun Liu, Marta Milo, Neil D. Lawrence, Magnus Rattray }
\seealso{Related methods \code{\link{pumaComb}}, \code{\link{mmgmos}} and \code{\link{pplr}}}
\examples{
  data(exampleE)
  data(exampleStd)
  r<-bcomb(exampleE,exampleStd,replicates=c(1,1,1,2,2,2),method="map")
}
\keyword{ manip }
\keyword{ models }
