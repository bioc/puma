\name{plotWhiskers}
\alias{plotWhiskers}
\title{Standard errors whiskers plot}
\description{
A plot showing error bars for genes of interest.
}
\usage{
plotWhiskers(
	eset
,	comparisons=c(1,2)
,	sortMethod = c("logRatio", "PPLR")
,	numGenes=50
,	xlim
,	main = "PUMA Whiskers plot"
,	highlightedGenes=NULL
)
}
\arguments{
	\item{eset}{ An object of class \code{\link[Biobase]{ExpressionSet}}. }
	\item{comparisons}{ A 2-element integer vector specifying the columns of data to be compared. }
	\item{sortMethod}{ The method used to sort the genes. "logRatio" is fold change. PPLR is Probability of Positive Log Ratio (as determined by the \code{\link{pumaDE}} method). }
	\item{numGenes}{ ~~Describe \code{numGenes} here~~ }
	\item{xlim}{ ~~Describe \code{xlim} here~~ }
	\item{main}{ ~~Describe \code{main} here~~ }
	\item{highlightedGenes}{ ~~Describe \code{highlightedGenes} here~~ }
}
\value{
	This function has no return value. The output is the plot created.
}
\author{ Richard D. Pearson }
\seealso{Related method \code{\link{pumaDE}}}
\keyword{hplot}
