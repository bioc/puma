\name{pumaCombImproved}
\alias{pumaCombImproved}
\title{Combining replicates for each condition with the true gene expression  }
\description{
This function calculates the combined (from replicates) signal for each condition using Bayesian models, which are added a hidden variable to represent the true expression for each gene on each chips. The inputs are gene expression levels and the probe-level standard deviations associated with expression measurements for each gene on each chip. The outputs include gene expression levels and standard deviation for each condition.
}
\usage{
pumaCombImproved(
	eset
,	design.matrix=NULL
,	numOfChunks=1000
,       maxOfIterations=200
,	save_r=FALSE
,	cl=NULL
,	parallelCompute=FALSE
)
}
\arguments{
	\item{eset}{ An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. }
	\item{design.matrix}{ A design matrix. }
	\item{numOfChunks}{ An integer defining how many chunks the data is divided into before processing. There is generally no need to change the default value. }
	\item{maxOfIterations}{ The maximum number of iterations controls the convergence. }
	\item{save_r}{ Will save an internal variable \code{r} to a file. Used for debugging purposes. }
\item{cl}{ A "cluster" object. See
	\code{\link[snow:snow-startstop]{makeCluster}} function from
	\pkg{snow} package for more details (if available).} 
	\item{parallelCompute}{ Boolean identifying whether processing in parallel should occur. }
}
\details{
It is generally recommended that data is normalised prior to using this function. Note that the default behaviour of mmgmos is to normalise data so this shouldn't generally be an issue. See the function \code{\link{pumaNormalize}} for more details on normalisation.

The maxOfIterations is used to control the maximum number of the iterations in the EM algorithm. You can change the number of maxOfIterations, but the best value of the maxOfIterations is from 200 to 1000, and should be set 200 at least. The default value is 200. 
}
\value{
  The result is an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object.
}
\references{Gelman,A., Carlin,J.B., Stern,H.S., Rubin,D.B., Bayesian data analysis. London: Chapman & Hall; 1995.

Zhang,L. and Liu,X. (2009) An improved probabilistic model for finding differential gene expression, technical report available request.
the 2nd BMEI 17-19 oct. 2009. Tianjin. China.

Liu,X., Milo,M., Lawrence,N.D. and Rattray,M. (2006)  Probe-level variances improve accuracy in detecting differential
gene expression, Bioinformatics, 22(17):2107-13.}
\author{ Li Zhang, Xuejun Liu }
\seealso{Related methods \code{\link{pumaNormalize}}, \code{\link{hcomb}}, \code{\link{mmgmos}} and \code{\link{pumaDE}}}
\examples{
	#	Next 4 lines commented out to save time in package checks, and saved version used
    # if (require(affydata)) {
	#	data(Dilution)
	#	eset_mmgmos <- mmgmos(Dilution)
	# }
	data(eset_mmgmos)

	#	Next line shows that eset_mmgmos has 4 arrays, each of which is a different
	#   condition (the experimental design is a 2x2 factorial, with both liver and
	#	scanner factors)
	pData(eset_mmgmos)
	
	#	Next line shows expression levels of first 3 probe sets
	exprs(eset_mmgmos)[1:3,]

	#	Next line used so eset_mmgmos only has information about the liver factor
	#	The scanner factor will thus be ignored, and the two arrays of each level
	#	of the liver factor will be treated as replicates
	pData(eset_mmgmos) <- pData(eset_mmgmos)[,1,drop=FALSE]

	#	To save time we'll just use 100 probe sets for the example
	eset_mmgmos_100 <- eset_mmgmos[1:100,]
	eset_combimproved <- pumaCombImproved(eset_mmgmos_100)
	
	#	We can see that the resulting ExpressionSet object has just two conditions
	#	and 1 expression level for each condition
	pData(eset_combimproved)
	exprs(eset_combimproved)[1:3,]
}
\keyword{ manip }
\keyword{ models }
