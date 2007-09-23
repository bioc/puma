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
,	parallelCompute=if(
		"Rmpi" \%in\% installed.packages() & "snow" \%in\% installed.packages()
	)
	as.logical(length(grep("origin",system("lamnodes",TRUE,TRUE))))
	else FALSE
)
}
\arguments{
	\item{eset}{ An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. }
	\item{design.matrix}{ A design matrix. }
	\item{method}{ Method "map" uses MAP of a hierarchical Bayesion model with Gamma prior on the between-replicate variance (Gelman et.al. p.285) and shares the same variance across conditions. This method is fast and suitable for the case where there are many conditions.

Method "em" uses variational inference of the same hierarchical Bayesian model as in method "map" but with conjugate prior on between-replicate variance and shares the variance across conditions. This is generaly much slower than "map", but is recommended where there are few conditions (as is usually the case). }
	\item{numOfChunks}{ An integer defining how many chunks the data is divided into before processing. There is generally no need to change the default value. }
	\item{save_r}{ Will save an internal variable \code{r} to a file. Used for debugging purposes. }
	\item{cl}{ A "cluster" object. See \code{\link[snow]{makeCluster}} function from \pkg{snow} package for more details (if available).}
	\item{parallelCompute}{ Boolean identifying whether processing in parallel should occur. }
}
\details{
It is generally recommended that data is normalised prior to using this function. See the function \code{\link{pumaNormalize}}.
}
\value{
  The result is an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object.
}
\references{Gelman,A., Carlin,J.B., Stern,H.S., Rubin,D.B., Bayesian data analysis. London: Chapman & Hall; 1995.

Liu,X., Milo,M., Lawrence,N.D. and Rattray,M. (2005) Probe-level variances improve accuracy in detecting differential
gene expression, technical report available upon request.}
\author{ Xuejun Liu, Marta Milo, Neil D. Lawrence, Magnus Rattray }
\seealso{Related methods \code{\link{pumaNormalize}}, \code{\link{bcomb}}, \code{\link{mmgmos}} and \code{\link{pumaDE}}}
\examples{
	data(affybatch.example)
	varLabels(affybatch.example) <- c("level")
	pData(affybatch.example) <- data.frame("level"=c("twenty","twenty","ten")
	    , "batch"=c("A","B","A"), row.names=rownames(pData(affybatch.example)))
	eset_mmgmos <- mmgmos(affybatch.example)
	pData(eset_mmgmos)
	exprs(eset_mmgmos)[1:3,]
	eset_comb <- pumaComb(eset_mmgmos)
	pData(eset_comb)
	exprs(eset_comb)[1:3,]
}
\keyword{ manip }
\keyword{ models }
