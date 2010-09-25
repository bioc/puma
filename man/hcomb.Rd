\name{hcomb}
\alias{hcomb}
\title{Combining replicates for each condition with the true gene expression }
\description{
 This function calculates the combined (from replicates) signal for each condition using Bayesian models, which are added a hidden variable to represent the true expression for each gene on each chips. The inputs are gene expression levels and the probe-level standard deviations associated with expression measurements for each gene on each chip. The outputs include gene expression levels and standard deviation for each condition.
}
\usage{
hcomb(e, se, replicates, max_num=c(200,500,1000),gsnorm=FALSE, eps=1.0e-6)
}
\arguments{
  \item{e}{ a data frame containing the expression level for each gene on each chip. }
  \item{se}{ a data frame containing the standard deviation of gene expression levels. }
  \item{replicates}{ a vector indicating which chip belongs to which condition. }  
  \item{max_num}{integer. The maximum number of iterations controls the convergence.  }
  \item{gsnorm}{ logical specifying whether do global scaling normalisation or not. }
  \item{eps}{ a numeric, optimisation parameter. }
}
\details{
  Each element in replicate represents the condition of the chip which is in the same column order as 
  in the expression and standard deviation matrix files.
  
  The max\_num is used to control the maximum number of the iterations in the EM algorithm. 
  The best value of the max\_num is from 200 to 1000, 
  and should be set 200 at least. The default value is 200.
   
  }
\value{
  The result is a data frame with components named 'M1', 'M2', and so on, which represent the 
  mean expression values for condition 1, condition 2, and so on. It also has components named 
  'Std1', 'Std2', and so on, which represent the standard deviation of the gene expression values 
  for condition 1, condtion 2, and so on.
}
\references{Gelman,A., Carlin,J.B., Stern,H.S., Rubin,D.B., Bayesian data analysis. London: Chapman & Hall; 1995.

Zhang,L. and Liu,X. (2009) An improved probabilistic model for finding differential gene expression, technical report available request.

Liu,X., Milo,M., Lawrence,N.D. and Rattray,M. (2006)  Probe-level variances improve accuracy in detecting differential
gene expression, Bioinformatics, 22(17):2107-13.
}
\author{ Li Zhang, Xuejun Liu }
\seealso{ Related method  \code{\link{pumaCombImproved}}, \code{\link{mmgmos}} and \code{\link{pplr}}}
\examples{
  data(exampleE)
  data(exampleStd)
  r<-hcomb(exampleE,exampleStd,replicates=c(1,1,1,2,2,2))
}
\keyword{ manip }
\keyword{ models }
