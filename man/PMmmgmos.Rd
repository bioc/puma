 \name{PMmmgmos}
\alias{PMmmgmos}
\title{ Multi-chip modified gamma Model for Oligonucleotide Signal using only PM probe intensities}
\description{
  This function converts an object of class \code{\link[oligoClasses:FeatureSet-class]{FeatureSet}} into an object of class
  \code{\link{exprReslt}} using the Multi-chip modified gamma Model for Oligonucleotide Signal
  (PMmulti-mgMOS). This method uses only PM probe intensites. This function obtains confidence of measures, standard deviation and 5,
  25, 50, 75 and 95 percentiles, as well as the estimated expression levels.
}
\usage{
PMmmgmos(
	object
,	background=TRUE
,	replaceZeroIntensities=TRUE
,	gsnorm=c("median", "none", "mean", "meanlog")
,	savepar=FALSE
,	eps=1.0e-6
,	addConstant = 0
)
}
\arguments{
  \item{object}{an object of \code{\link[oligoClasses:FeatureSet-class]{FeatureSet}}}
  \item{background}{Logical value. If \code{TRUE}, perform background correction before applying PMmmgmos.}
  \item{replaceZeroIntensities}{Logical value. If \code{TRUE}, replace 0 intensities with 1 before applying PMmmgmos.}
  \item{gsnorm}{character. specifying the algorithm of global scaling normalisation.}
  \item{savepar}{Logical value. If \code{TRUE} the estimated parameters of the model are saved in file par\_pmmmgmos.txt  }
  \item{eps}{Optimisation termination criteria.}
  \item{addConstant}{numeric. This is an experimental feature and should not generally be changed from the default value. }
}
\details{
  The obtained expression measures are in log base 2 scale.

  The algorithms of global scaling normalisation can be one of "median", "none", "mean", "meanlog".
  "mean" and "meanlog" are mean-centered normalisation on raw scale and log scale respectively, and "median"
  is median-centered normalisation. "none" will result in no global scaling normalisation being applied.

  There are n+2 columns in file par\_pmmmgmos.txt, n is the number of chips. The first n columns are 'alpha' values for n chips,
  column n+1 is 'c' values and the final column is values for 'd'.
}
\value{
  An object of class \code{exprReslt}.
}
\references{Liu,X., Milo,M., Lawrence,N.D. and Rattray,M. (2005) A tractable probabilistic model for
Affymetrix probe-level analysis across multiple chips, Bioinformatics 21: 3637-3644.

Milo,M., Niranjan,M., Holley,M.C., Rattray,M. and Lawrence,N.D. (2004) A probabilistic approach for summarising
oligonucleotide gene expression data, technical report available upon request.

Milo,M., Fazeli,A., Niranjan,M. and Lawrence,N.D. (2003) A probabilistic model for the extractioin of expression
 levels from oligonucleotide arrays, Biochemical Society Transactions, 31: 1510-1512.

Peter Spellucci. DONLP2 code and accompanying documentation. Electronically available via
http://plato.la.asu.edu/donlp2.html}

\author{Xuejun Liu, Zhenzhu Gao, Magnus Rattray, Marta Milo, Neil D. Lawrence}

\seealso{ Related class \code{\link{exprReslt-class}} and related method \code{\link{mgmos}}}
\examples{
## Code commented out to speed up checks
## load example data from package pumadata
 #if (require(pumadata)&&require(puma)){
   # data(oligo.estrogen)
## use method PMmmgMOS to calculate the expression levels and related confidence
##of the measures for the example data
 #   eset<-PMmmgmos(oligo.estrogen,gsnorm="none")
#}
}
\keyword{manip}
