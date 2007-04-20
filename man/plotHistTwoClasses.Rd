\name{plotHistTwoClasses}
\alias{plotHistTwoClasses}
\title{Stacked histogram plot of two different classes}
\description{
Stacked histogram plot of two different classes
}
\usage{
plotHistTwoClasses(
	scores
,	class1Elements
,	class2Elements
,	space=0
,	col=c("white", "grey40")
,	xlab="PPLR"
,	ylab="Number of genes"
,	ylim=NULL
,	las=0 # axis labels all perpendicular to axes
,	legend=c("non-spike-in genes", "spike-in genes")
,	inset=0.05
,	minScore=0
,	maxScore=1
,	numOfBars=20
,	main=NULL
)
}
\arguments{
	\item{scores}{ A numeric vector of scores (e.g. from the output of \code{\link{pumaDE}}) }
	\item{class1Elements}{ Boolean vector, TRUE if element is in first class }
	\item{class2Elements}{ Boolean vector, TRUE if element is in second class }
	\item{space}{ Numeric. x-axis distance between bars }
	\item{col}{ Colours for the two different classes }
	\item{xlab}{ Title for the x-axis }
	\item{ylab}{ Title for the y-axis }
	\item{ylim}{ 2-element numeric vector showing minimum and maximum values for y-axis. }
	\item{las}{ See \code{\link[graphics]{par}}. Default of 0 means axis labels all perpendicular to axes. }
	\item{legend}{ 2-element string vector giving text to appear in legend for the two classes. }
	\item{inset}{ See \code{\link[graphics]{legend}} }
	\item{minScore}{ Numeric. Minimum score to plot. }
	\item{maxScore}{ Numeric. Maximum score to plot. }
	\item{numOfBars}{ Integer. Number of bars to plot. }
	\item{main}{ String. Main title for the plot. }
}
\value{
	This function has no return value. The output is the plot created.
}
\author{ Richard D. Pearson }
\examples{
	class1 <- rnorm(1000,0.2,0.1)
	class2 <- rnorm(1000,0.6,0.2)
	class1[which(class1<0)] <- 0
	class1[which(class1>1)] <- 1
	class2[which(class2<0)] <- 0
	class2[which(class2>1)] <- 1
	scores <- c(class1, class2)
	class1elts <- c(rep(TRUE,1000), rep(FALSE,1000))
	class2elts <- c(rep(FALSE,1000), rep(TRUE,1000))
	plotHistTwoClasses(scores, class1elts, class2elts, ylim=c(0,300))
}
\keyword{hplot}
