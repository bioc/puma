\name{pumaFull}
\alias{pumaFull}
\title{Perform a full PUMA analysis}
\description{
Full analysis including pumaPCA and mmgmos/pumaDE vs rma/limma comparison
}
\usage{
pumaFull (
	ExpressionFeatureSet = NULL
,	data_dir = getwd()
,	load_ExpressionFeatureSet = FALSE
,	calculate_eset = TRUE
,	calculate_pumaPCAs = TRUE
,	calculate_bcomb = TRUE
,	mmgmosComparisons = FALSE
)

}
\arguments{
	\item{ExpressionFeatureSet}{ An object of class \code{\link[oligoClasses:FeatureSet-class]{FeatureSet}}. }
	\item{data_dir}{ A character string specifying where data files are stored. }
	\item{load_ExpressionFeatureSet}{ Boolean. Load a pre-existing ExpressionFeatureSet object? Note that this has to be named "ExpressionFeatureSet.rda" and be in the \code{data_dir} directory. }
	\item{calculate_eset}{ Boolean. Calculate ExpressionSet from \code{ExpressionFeatureSet} object? If FALSE, files named "eset\_mmgmos.rda" and "eset\_rma.rda" must be available in the \code{data\_dir} directory. }
	\item{calculate_pumaPCAs}{ Boolean. Calculate pumaPCA from \code{eset\_mmgmos} object? If FALSE, a file named "pumaPCA\_results.rda" must be available in the \code{data\_dir} directory. }
	\item{calculate_bcomb}{ Boolean. Calculate pumaComb from \code{eset\_mmgmos} object? If FALSE, files named "eset\_comb.rda" and "eset\_normd\_comb.rda" must be available in the \code{data\_dir} directory. }
	\item{mmgmosComparisons}{ Boolean. If TRUE, will compare mmgmos with default settings, with mmgmos used with background correction. }
}
\value{
  No return values. Various objects are saved as .rda files during the execution of this function, and various PDF files are created.
}
\author{ Richard D. Pearson }
\seealso{Related methods \code{\link{pumaDE}}, \code{\link{createDesignMatrix}} and \code{\link{createContrastMatrix}}}
\examples{
##	Code commented out to ensure checks run quickly
#	if (require(pumadata)) data(oligo.estrogen)
#	pumaFull(oligo.estrogen)
}
\keyword{manip}
