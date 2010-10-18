\name{pumaComb}
\alias{pumaComb}
\title{Combining replicates for each condition}
\description{
This function calculates the combined (from replicates) signal for each condition using Bayesian models. The inputs are gene expression levels and the probe-level standard deviations associated with expression measurements for each gene on each chip. The outputs include gene expression levels and standard deviation for each condition.
}
\usage{
pumaComb(
	eset
,	design.matrix=NULL
,	method="em"
,	numOfChunks=1000
,	save_r=FALSE
,	cl=NULL
,	parallelCompute=FALSE
)
}
\arguments{
	\item{eset}{ An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. }
	\item{design.matrix}{ A design matrix. }
	\item{method}{ Method "map" uses MAP of a hierarchical Bayesion model with Gamma prior on the between-replicate variance (Gelman et.al. p.285) and shares the same variance across conditions. This method is fast and suitable for the case where there are many conditions.

Method "em" uses variational inference of the same hierarchical Bayesian model as in method "map" but with conjugate prior on between-replicate variance and shares the variance across conditions. This is generaly much slower than "map", but is recommended where there are few conditions (as is usually the case). }
	\item{numOfChunks}{ An integer defining how many chunks the data is divided into before processing. There is generally no need to change the default value. }
	\item{save_r}{ Will save an internal variable \code{r} to a file. Used for debugging purposes. }
\item{cl}{ A "cluster" object. See
	\code{\link[snow:snow-startstop]{makeCluster}} function from
	\pkg{snow} package for more details (if available).} 
	\item{parallelCompute}{ Boolean identifying whether processing in parallel should occur. }
}
\details{
It is generally recommended that data is normalised prior to using this function. Note that the default behaviour of mmgmos is to normalise data so this shouldn't generally be an issue. See the function \code{\link{pumaNormalize}} for more details on normalisation.
}
\value{
  The result is an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object.
}
\references{Gelman,A., Carlin,J.B., Stern,H.S., Rubin,D.B., Bayesian data analysis. London: Chapman & Hall; 1995.

Liu,X., Milo,M., Lawrence,N.D. and Rattray,M. (2006) Probe-level variances improve accuracy in detecting differential
gene expression, Bioinformatics, 22:2107-2113.
}
\author{ Xuejun Liu, Marta Milo, Neil D. Lawrence, Magnus Rattray }
\seealso{Related methods \code{\link{pumaNormalize}}, \code{\link{bcomb}}, \code{\link{mmgmos}} and \code{\link{pumaDE}}}
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
	eset_comb <- pumaComb(eset_mmgmos_100)
	
	#	We can see that the resulting ExpressionSet object has just two conditions
	#	and 1 expression level for each condition
	pData(eset_comb)
	exprs(eset_comb)[1:3,]
}
\keyword{ manip }
\keyword{ models }
