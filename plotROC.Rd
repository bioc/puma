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
,	includedProbesets=1:length(truthValues)
,	legendTitles=1:length(scoresList)
,	main = "PUMA ROC plot"
,	lty = 1:length(scoresList)
,	col = rep(1,length(scoresList))
,	lwd = rep(1,length(scoresList))
,	yaxisStat = "tpr"
,	xaxisStat = "fpr"
,	downsampling = 100
,	showLegend = TRUE
,	showAUC = TRUE
,	...
)
}
\arguments{
	\item{scoresList}{ A list, each element of which is a numeric vector of scores. }
	\item{truthValues}{ A boolean vector indicating which scores are True Positives. }
	\item{includedProbesets}{ A vector of indices indicating which scores (and truthValues) are to be used in the calculation. The default is to use all, but a subset can be used if, for example, you only want a subset of the probesets which are not True Positives to be treated as False Positives. }
	\item{legendTitles}{ Vector of names to appear in legend. }
	\item{main}{ Main plot title }
	\item{lty}{ Line types. }
	\item{col}{ Colours. }
	\item{lwd}{ Line widths. }
	\item{yaxisStat}{Character string identifying what is to be plotted on the y-axis. The default is "tpr" for True Positive Rate. See \code{\link[ROCR:performance]{performance}} function from \pkg{ROCR} package.}
	\item{xaxisStat}{Character string identifying what is to be plotted on the x-axis. The default is "fpr" for False Positive Rate. See \code{\link[ROCR:performance]{performance}} function from \pkg{ROCR} package.}
\item{downsampling}{See details for
	  \code{\link[ROCR:plot-methods]{plot.performance}} from the
	  \pkg{ROCR} package. } 
	\item{showLegend}{Boolean. Should legend be displayed?}
	\item{showAUC}{Boolean. Should AUC values be included in legend?}
	\item{\dots}{ Other parameters to be passed to \code{plot}. }
}
\value{
	This function has no return value. The output is the plot created.
}
\author{ Richard D. Pearson }
\seealso{Related method \code{\link{calcAUC}}}
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
