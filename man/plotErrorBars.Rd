\name{plotErrorBars}
\alias{plotErrorBars}
\title{Plot mean expression levels and error bars for one or more probesets}
\description{
This produces plots of probesets of interest.
}
\usage{
plotErrorBars(
	eset
,	probesets = if(dim(exprs(eset))[1] <= 12) 1:dim(exprs(eset))[1] else 1
,	arrays = 1:dim(pData(eset))[1] # default is to use all
,	xlab = paste(colnames(pData(eset))[1:numOfFactorsToUse(eset)], collapse=":")
,	ylab = "Expression Estimate"
,	xLabels = apply(
                  as.matrix(pData(eset)[arrays,1:numOfFactorsToUse(eset)])
                , 1
                , function(mat){paste(mat, collapse=":")}
                )
,	ylim = NA
,	numOfSEs = qnorm(0.975)
,	globalYlim = FALSE # Not yet implemented!
,	plot_cols = NA
,	plot_rows = NA
,	featureNames = NA
,	showGeneNames = FALSE
,	showErrorBars = if(
					length(assayDataElement(eset,"se.exprs"))==0 ||
					length(assayDataElement(eset,"se.exprs")) == sum(is.na(assayDataElement(eset,"se.exprs")))
					) FALSE else TRUE
,	plotColours = FALSE
,	log.it = if(max(exprs(eset)) > 32) TRUE else FALSE
,	eset_comb = NULL
,	jitterWidth = NA
,	qtpcrData = NULL
, ...
)
}
\arguments{
	\item{eset}{ An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. This is the main object being plotted. }
	\item{probesets}{ A vector of integers indicating the probesets to be plotted. These integers refer to the row numbers of the \code{eset}. }
	\item{arrays}{ A vector of integers indicating the arrays to be shown on plots. }
	\item{xlab}{ Character string of title to appear on x-axis }
	\item{ylab}{ Character string of title to appear on y-axis }
	\item{xLabels}{ Vector of strings for labels of individual points on x-axis. }
	\item{ylim}{ 2-element numeric vector showing minimum and maximum values for y-axis. }
	\item{numOfSEs}{ Numeric indicating the scaling for the error bars. The default value give error bars that include 95\% of expected values. }
	\item{globalYlim}{ Not yet implemented! }
	\item{plot_cols}{ Integer specifying number of columns for multi-figure plot. }
	\item{plot_rows}{ Integer specifying number of rows for multi-figure plot. }
	\item{featureNames}{ A vector of strings for \code{featureNames} (Affy IDs). This is an alternative (to the \code{probesets} argument) way of specifying probe sets. }
	\item{showGeneNames}{ Boolean indicating whether to use Affy IDs as titles for each plot. }
	\item{showErrorBars}{ Boolean indicating whether error bars should be shown on plots. }
	\item{plotColours}{ A vector of colours to plot. }
	\item{log.it}{ Boolean indicating whether expression values should be logged. }
	\item{eset_comb}{ An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. This is a secondary object to be plotted on the same charts as \code{eset}. This should be an object created using \code{\link{pumaComb}} and \code{\link{pumaCombImproved}} which holds the values created by combining information from the replicates of each condition. }
	\item{jitterWidth}{ Numeric indicating the x-axis distance between replicates of the same condition. }
	\item{qtpcrData}{ A 2-column matrix of qRT-PCR values (or other data to be plotted on the same charts). }
	\item{\dots}{ Additional arguments to be passed to \code{\link[graphics]{plot}}. }
}
\value{
	This function has no return value. The output is the plot created.
}
\author{ Richard D. Pearson }
\examples{
	#	Next 4 lines commented out to save time in package checks, and saved version used
    # if (require(affydata)) {
	#	data(Dilution)
	#	eset_mmgmos <- mmgmos(Dilution)
	# }
	data(eset_mmgmos)
	plotErrorBars(eset_mmgmos)
	plotErrorBars(eset_mmgmos,1:6)
}

\keyword{hplot}
