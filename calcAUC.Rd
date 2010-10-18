\name{calcAUC}
\alias{calcAUC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate Area Under Curve (AUC) for a standard ROC plot. }
\description{
	Calculates the AUC values for one or more ROC plots.
}
\usage{
calcAUC(scores, truthValues, includedProbesets = 1:length(truthValues))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scores}{ A vector of scores. This could be, e.g. one of the columns of the statistics of a \code{\link{DEResult}} object. }
  \item{truthValues}{ A boolean vector indicating which scores are True Positives. }
  \item{includedProbesets}{ A vector of indices indicating which scores (and truthValues) are to be used in the calculation. The default is to use all, but a subset can be used if, for example, you only want a subset of the probesets which are not True Positives to be treated as False Positives. }
}
\value{
  A single number which is the AUC value.
}
\author{ Richard D. Pearson }
\seealso{Related methods \code{\link{plotROC}} and \code{\link{numFP}}.}
\examples{
class1a <- rnorm(1000,0.2,0.1)
class2a <- rnorm(1000,0.6,0.2)
class1b <- rnorm(1000,0.3,0.1)
class2b <- rnorm(1000,0.5,0.2)
scores_a <- c(class1a, class2a)
scores_b <- c(class1b, class2b)
classElts <- c(rep(FALSE,1000), rep(TRUE,1000))
print(calcAUC(scores_a, classElts))
print(calcAUC(scores_b, classElts))

}
\keyword{manip}
