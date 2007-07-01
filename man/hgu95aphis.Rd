\name{hgu95aphis}
\alias{hgu95aphis}
\docType{data}
\title{ Estimated parameters of the distribution of phi}
\description{
  The pre-estimated parameters of log normal distribution of \eqn{\phi}, 
  which is the fraction of specific signal binding to mismatch probe.
}
\usage{data(hgu95aphis)}
\format{
  The format is:
 num [1:3]  0.171 -1.341  0.653
}
\details{
  The current values of hgu95aphis are estimated from Affymetrix spike-in data sets.
  It was loaded in the method \code{"mmgmos"}.
  
  hgu95aphis[1:3] is respectively the mode, mean and variance of the log normal distribution 
  of \eqn{\phi}, and hgu95aphis[1] is also the intial value of \eqn{\phi} in the model optimisation.
}
%\source{}
%\references{}
\keyword{datasets}
