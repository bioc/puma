\name{numFP}
\alias{numFP}
\title{ Number of False Positives for a given proportion of True Positives. }
\description{
  Often when evaluating a differential expression method, we are interested in how well a classifier performs for very small numbers of false positives. This method gives one way of calculating this, by determining the number of false positives for a set proportion of true positives.
}
\usage{
numFP(scores, truthValues, TPRate = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scores}{ A vector of scores. This could be, e.g. one of the columns of the statistics of a \code{\link{DEResult}} object. }
  \item{truthValues}{ A boolean vector indicating which scores are True Positives. }
  \item{TPRate}{ A number between 0 and 1 identify the proportion of true positives for which we wish to determine the number of false positives. }
}
\value{
  An integer giving the number of false positives.
}
\author{ Richard D. Pearson }
\seealso{Related methods \code{\link{plotROC}} and \code{\link{calcAUC}}.}
\examples{
class1a <- rnorm(1000,0.2,0.1)
class2a <- rnorm(1000,0.6,0.2)
class1b <- rnorm(1000,0.3,0.1)
class2b <- rnorm(1000,0.5,0.2)
scores_a <- c(class1a, class2a)
scores_b <- c(class1b, class2b)
classElts <- c(rep(FALSE,1000), rep(TRUE,1000))
print(numFP(scores_a, classElts))
print(numFP(scores_b, classElts))
}
\keyword{manip}
