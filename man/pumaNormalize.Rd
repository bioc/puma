\name{pumaNormalize}
\alias{pumaNormalize}
\title{Normalize an ExpressionSet}
\description{
This is used to apply a scaling normalization to set of arrays. This normalization can be at the array scale (thus giving all arrays the same mean or median), or at the probeset scale (thus giving all probesets the same mean or median).

It is generally recommended that the default option (median array scaling) is used after running \code{\link{mmgmos}} and before running \code{\link{pumaComb}} and/or \code{\link{pumaDE}}. There are however, situations where this might not be the recommended, for example in time series experiments where it is expected than there will be general up-regulation or down-regulation in overall gene expression levels between time points.
}
\usage{
pumaNormalize(
	eset
,	arrayScale = c("median", "none", "mean", "meanlog")
,	probesetScale = c("none", "mean", "median")
,	probesetNormalisation = NULL
,	replicates = list(1:dim(exprs(eset))[2])
)
}
\arguments{
	\item{eset}{ An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. }
	\item{arrayScale}{ A method of scale normalisation at the array level. }
	\item{probesetScale}{ A method of scale normalisation at the probe set level. }
	\item{probesetNormalisation}{ If not NULL normalises the expression levels to have zero mean and adjusts the variance of the gene expression according to the zero-centered normalisation. }
	\item{replicates}{ List of integer vectors indicating which arrays are replicates. }
}
\value{
An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} holding the normalised data.
}
\author{ Richard D. Pearson }
\seealso{Methods \code{\link{mmgmos}}, \code{\link{pumaComb}} and \code{\link{pumaDE}}}
\examples{
	#	Next 4 lines commented out to save time in package checks, and saved version used
    # if (require(affydata)) {
	#	data(Dilution)
	#	eset_mmgmos <- mmgmos(Dilution)
	# }
	data(eset_mmgmos)
	apply(exprs(eset_mmgmos),2,median)
	eset_mmgmos_normd <- pumaNormalize(eset_mmgmos)
	apply(exprs(eset_mmgmos_normd),2,median)
}

\keyword{manip}
