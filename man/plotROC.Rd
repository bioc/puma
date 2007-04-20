\name{plotROC}
\alias{plotROC}
\title{Receiver Operator Characteristic (ROC) plot}
\description{
Plots a Receiver Operator Characteristic (ROC) curve.
}
\usage{
plotROC(
	scoresList
,	truthValues
,	legendTitles=1:length(scoresList)
,	main = "PUMA ROC plot"
,	lty = 1:length(scoresList)
,	col = rep(1,length(scoresList))
,	...
)
}
\arguments{
	\item{scoresList}{ A list, each element of which is a numeric vector of scores }
	\item{truthValues}{ A boolean vector indicating which scores are True Positives }
	\item{legendTitles}{ Vector of names to appear in legend. }
	\item{main}{ Main plot title }
	\item{lty}{ Line types. }
	\item{col}{ Colours. }
	\item{\dots}{ Other parameters to be passed to \code{plot}. }
}
\value{
	This function has no return value. The output is the plot created.
}
\author{ Richard D. Pearson }
\examples{
	class1a <- rnorm(1000,0.2,0.1)
	class2a <- rnorm(1000,0.6,0.2)
	class1b <- rnorm(1000,0.3,0.1)
	class2b <- rnorm(1000,0.5,0.2)
	scores_a <- c(class1a, class2a)
	scores_b <- c(class1b, class2b)
	scores <- list(scores_a, scores_b)
	classElts <- c(rep(FALSE,1000), rep(TRUE,1000))
	plotROC(scores, classElts)
}
\keyword{hplot}
