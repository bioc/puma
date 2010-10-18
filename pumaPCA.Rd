\name{pumaPCA}
\alias{pumaPCA}
\alias{pumaPCAEstep}
\alias{pumaPCALikelihoodBound}
\alias{pumaPCALikelihoodCheck}
\alias{pumaPCANewtonUpdateLogSigma}
\alias{newtonStep}
\alias{pumaPCARemoveRedundancy}
\alias{pumaPCASigmaGradient}
\alias{pumaPCASigmaObjective}
\alias{pumaPCAUpdateCinv}
\alias{pumaPCAUpdateM}
\alias{pumaPCAUpdateMu}
\alias{pumaPCAUpdateW}
\title{PUMA Principal Components Analysis}
\description{
This function carries out principal components analysis (PCA), taking into
account not only the expression levels of genes, but also the variability
in these expression levels.

The various other pumaPCA... functions are called during the execution of pumaPCA
}
\usage{
pumaPCA(
    eset
,   latentDim           = 	if(dim(exprs(eset))[2] <= 3)
 								dim(exprs(eset))[[2]]-1
							else
								3
,   sampleSize          =	if(dim(exprs(eset))[1] <= 1000)
								 dim(exprs(eset))[[1]]
							else
								1000	## Set to integer or FALSE for all
,   initPCA             =	TRUE	## Initialise parameters with PCA
,   randomOrder         =	FALSE	## Update parameters in random order
,   optimMethod         =	"BFGS"  ## ?optim for details of methods
,   stoppingCriterion   =	"deltaW"## can also be "deltaL"
,   tol                 =	1e-3	## Stop when delta update < this
,   stepChecks          =	FALSE	## Check likelihood after each update?
,   iterationNumbers    =	TRUE	## Show iteration numbers?
,   showUpdates         =	FALSE	## Show values after each update?
,   showTimings         =	FALSE	## Show timings after each update?
,   showPlot            =	FALSE	## Show projection plot after each update?
,   maxIters            =	500		## Number of EM iterations.
,   transposeData       =	FALSE	## Transpose eset matrices?
,   returnExpectations  =	FALSE
,   returnData          =	FALSE
,   returnFeedback      =	FALSE
,	pumaNormalize		=	TRUE
)
}
\arguments{
	\item{eset}{ An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. }
	\item{latentDim}{ An integer specifying the number of latent dimensions (kind of like the number of principal components). }
	\item{sampleSize}{ An integer specifying the number of probesets to sample (default is 1000), or FALSE, meaning use all the data. }
	\item{initPCA}{ A boolean indicating whether to initialise using standard PCA (the default, and generally quicker and recommended).  }
	\item{randomOrder}{ A boolean indicating whether the parameters should be updated in a random order (this is generally not recommended, and the default is FALSE). }
	\item{optimMethod}{ See ?optim for details of methods. }
	\item{stoppingCriterion}{ If set to "deltaW" will stop when W changes by less than \code{tol}. If "deltaL" will stop when L (lambda) changes by less than \code{tol}. }
	\item{tol}{ Tolerance value for \code{stoppingCriterion}. }
	\item{stepChecks}{ Boolean. Check likelihood after each update? }
	\item{iterationNumbers}{ Boolean. Show iteration numbers? }
	\item{showUpdates}{ Boolean. Show values after each update? }
	\item{showTimings}{ Boolean. Show timings after each update? }
	\item{showPlot}{ Boolean. Show projection plot after each update? }
	\item{maxIters}{ Integer. Maximum number of EM iterations. }
	\item{transposeData}{ Boolean. Transpose eset matrices? }
	\item{returnExpectations}{ Boolean. Return expectation values? }
	\item{returnData}{ Boolean. Return expectation data? }
	\item{returnFeedback}{ Boolean. Return feedback on progress of optimisation? }
	\item{pumaNormalize}{ Boolean. Normalise data prior to running algorithm (recommended)? }
}
\value{
  An object of class \code{\link{pumaPCARes}}
}
\author{ Richard D. Pearson }
\seealso{Related methods \code{\link{pumaDE}}, \code{\link{createDesignMatrix}} and \code{\link{createContrastMatrix}}}
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
\keyword{multivariate}
