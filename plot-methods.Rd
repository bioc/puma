\name{plot-methods}
\docType{methods}
\alias{plot.pumaPCARes}
\alias{plot-methods}
\alias{plot,pumaPCARes-method}
\alias{plot,pumaPCARes,missing-method}
\title{Plot method for pumaPCARes objects}
\description{
This is the method to plot objects of class pumaPCARes. It will produce a scatter plot of two of the principal components
}
\usage{
\S4method{plot}{pumaPCARes,missing}(..., firstComponent = 1, secondComponent = 2, useFilenames = FALSE, phenotype = pData(pumaPCARes@phenoData), legend1pos = "topright", legend2pos = "bottomright")
}

\arguments{
  \item{\dots}{Optional graphical parameters to adjust different components of the plot}
  \item{firstComponent}{Integer identifying which principal component to plot on the x-axis}
  \item{secondComponent}{Integer identifying which principal component to plot on the x-axis}
  \item{useFilenames}{Boolean. If TRUE then use filenames as plot points. Otherwise just use points.}
  \item{phenotype}{Phenotype information}
  \item{legend1pos}{String indicating where to put legend for first factor}
  \item{legend2pos}{String indicating where to put legend for second factor}
}
\examples{
	#	Next 4 lines commented out to save time in package checks, and saved version used
    # if (require(affydata)) {
	#	data(Dilution)
	#	eset_mmgmos <- mmgmos(Dilution)
	# }
	data(eset_mmgmos)
	pumapca_mmgmos <- pumaPCA(eset_mmgmos)
	plot(pumapca_mmgmos)
}
\keyword{methods}
\keyword{hplot}
