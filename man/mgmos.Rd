\name{mgmos}
\alias{mgmos}
\title{ modified gamma Model for Oligonucleotide Signal }
\description{
  This function converts an object of class \code{\link[oligoClasses:FeatureSet-class]{FeatureSet}} into an object of class 
  \code{\link{exprReslt}} using the modified gamma Model for Oligonucleotide Signal 
  (multi-mgMOS). This function obtains confidence of measures, standard deviation and 5, 
  25, 50, 75 and 95 percentiles, as well as the estimated expression levels.
}
\usage{
mgmos(
	object
,	background=FALSE
,	replaceZeroIntensities=TRUE
,	gsnorm=c("median", "none", "mean", "meanlog")
,	savepar=FALSE
,	eps=1.0e-6
)
}
\arguments{
	\item{object}{an object of \code{\link[oligoClasses:FeatureSet-class]{FeatureSet}}}
	\item{background}{Logical value. If \code{TRUE}, perform background correction before applying mmgmos.}
	\item{replaceZeroIntensities}{Logical value. If \code{TRUE}, replace 0 intensities with 1 before applying mmgmos.}
	\item{gsnorm}{character. specifying the algorithm of global scaling normalisation.}
	\item{savepar}{Logical value. If \code{TRUE} the estimated parameters of the model are saved in file par\_mmgmos.txt and phi\_mmgmos.txt. }
	\item{eps}{Optimisation termination criteria.}
}
\details{
  The obtained expression measures are in log base 2 scale.
  
  The algorithms of global scaling normalisation can be one of "median", "none", "mean", "meanlog".
  "mean" and "meanlog" are mean-centered normalisation on raw scale and log scale respectively, and "median" 
  is median-centered normalisation. "none" will result in no global scaling normalisation being applied.
  
  There are 4*n columns in file par\_mgmos.txt, n is the number of chips. Every 4 columns are parameters for a chip. 
  Among every 4 columns, the first one is for 'alpha' values, the 2nd one is for 'a' values, 
  The 3rd column is for 'c' and the final column is values for 'd'. 
}
\value{
  An object of class \code{exprReslt}.
}
\references{Liu,X., Milo,M., Lawrence,N.D. and Rattray,M. (2005) A tractable probabilistic model for 
Affymetrix probe-level analysis across multiple chips, Bioinformatics, 21:3637-3644.

Milo,M., Niranjan,M., Holley,M.C., Rattray,M. and Lawrence,N.D. (2004) A probabilistic approach for summarising 
oligonucleotide gene expression data, technical report available upon request.

Milo,M., Fazeli,A., Niranjan,M. and Lawrence,N.D. (2003) A probabilistic model for the extractioin of expression
 levels from oligonucleotide arrays, Biochemical Society Transactions, 31: 1510-1512.

Peter Spellucci. DONLP2 code and accompanying documentation. Electronically available via 
http://plato.la.asu.edu/donlp2.html}

\author{Xuejun Liu, Magnus Rattray, Marta Milo, Neil D. Lawrence}

\seealso{ Related class \code{\link{exprReslt-class}} and related method \code{\link{mmgmos}}}
\examples{
## Code commented out to speed up checks
## load example data from package affydata
# if (require(pumadata)&&require(puma)){ 
  #  data(oligo.estrogen)
# use method mgMOS to calculate the expression levels and related confidence 
# of the measures for the example data
  #  eset<-mgmos(oligo.estrogen,gsnorm="none")
#}
}
\keyword{manip}
