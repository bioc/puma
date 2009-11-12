\name{pplr}
\alias{pplr}
\title{ Probability of positive log-ratio }
\description{
WARNING - this function is generally not expected to be used, but is intended as an internal function. It is included for backwards compatibility with the \pkg{pplr} package, but may be deprecated and then hidden in future. Users should generally use \code{\link{pumaDE}} instead.

  This function calculates the probability of positive log-ratio (PPLR) between any two specified conditions in the input data, mean and standard deviation of gene expression level for each condition.
}
\usage{
pplr(e, control, experiment, sorted=TRUE)
}
\arguments{
  \item{e}{ a data frame containing the mean and standard deviation of gene expression 
  levels for each condition. }
  \item{control}{ an integer denoting the control condition. }
  \item{experiment}{ an integer denoting the experiment condition. }
  \item{sorted}{ Boolean. Should PPLR values be sorted by value? If FALSE, PPLR values are returned in same order as supplied.}
}
\details{
  The input of 'e' should be a data frame comprising of 2*n components, where n is 
  the number of conditions. The first 1,2,...,n components include the mean of gene expression 
  values for conditions 1,2,...,n, and the n+1, n+2,...,2*n components contain the standard deviation of 
  expression levels for condition 1,2,...,n.
}
\value{
  The return is a data frame. The description of the components are below.
  \item{index }{The original row number of genes.}
  \item{cM }{The mean expression levels under control condition.}
  \item{sM}{The mean expression levels under experiment condition.}
  \item{cStd}{The standard deviation of gene expression levels under control condition.}
  \item{sStd}{The standard deviation of gene expression levels under experiment condition.}
  \item{LRM}{The mean log-ratio between control and experiment genes.}
  \item{LRStd}{The standard deviation of log-ratio between control and experiment genes.}
  \item{stat}{A statistic value which is -mean/(sqrt(2)*standard deviation).}
  \item{PPLR}{Probability of positive log-ratio.}
}
\references{ Liu,X., Milo,M., Lawrence,N.D. and Rattray,M. (2006)  Probe-level variances improve accuracy in detecting differential
gene expression, Bioinformatics, 22(17):2107-13.}
\author{ Xuejun Liu, Marta Milo, Neil D. Lawrence, Magnus Rattray }
\seealso{ Related methods \code{\link{pumaDE}}, \code{\link{bcomb}} and \code{\link{hcomb}} }
\examples{
  data(exampleE)
  data(exampleStd)
  r<-bcomb(exampleE,exampleStd,replicates=c(1,1,1,2,2,2),method="map")
  p<-pplr(r,1,2)
}
\keyword{ manip }
\keyword{ models }
 
