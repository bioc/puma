\name{numTP}
\alias{numTP}
\title{ Number of True Positives for a given proportion of False Positives. }
\description{
  Often when evaluating a differential expression method, we are interested in how well a classifier performs for very small numbers of true positives. This method gives one way of calculating this, by determining the number of true positives for a set proportion of false positives.
}
\usage{
numTP(scores, truthValues, FPRate = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scores}{ A vector of scores. This could be, e.g. one of the columns of the statistics of a \code{\link{DEResult}} object. }
  \item{truthValues}{ A boolean vector indicating which scores are True Positives. }
  \item{FPRate}{ A number between 0 and 1 identify the proportion of flase positives for which we wish to determine the number of true positives. }
}
\value{
  An integer giving the number of true positives.
}
\author{ Richard D. Pearson }
\seealso{Related methods \code{\link{numFP}}, \code{\link{plotROC}} and \code{\link{calcAUC}}.}
\examples{
class1a <- rnorm(1000,0.2,0.1)
class2a <- rnorm(1000,0.6,0.2)
class1b <- rnorm(1000,0.3,0.1)
class2b <- rnorm(1000,0.5,0.2)
scores_a <- c(class1a, class2a)
scores_b <- c(class1b, class2b)
classElts <- c(rep(FALSE,1000), rep(TRUE,1000))
print(numTP(scores_a, classElts))
print(numTP(scores_b, classElts))
}
\keyword{manip}
