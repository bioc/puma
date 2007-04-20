\name{mmgmos}
\alias{mmgmos}
\title{ Multi-chip modified gamma Model for Oligonucleotide Signal }
\description{
  This function converts an object of class \code{AffyBatch} into an object of class 
  \code{exprReslt} using the Multi-chip modified gamma Model for Oligonucleotide Signal 
  (multi-mgMOS). This function obtains confidence of measures, standard deviation and 5, 
  25, 50, 75 and 95 percentiles, as well as the estimated expression levels.
}
\usage{
mmgmos(
	object
,	background=FALSE
,	replaceZeroIntensities=TRUE
,	gsnorm=c("none","mean","median","meanlog")
,	savepar=FALSE
,	eps=1.0e-6
,	orig.phis = FALSE
)
}
\arguments{
  \item{object}{an object of \code{\link[affy:AffyBatch-class]{AffyBatch}}}
  \item{background}{Logical value. If \code{TRUE}, perform background correction before applying mmgmos.}
  \item{replaceZeroIntensities}{Logical value. If \code{TRUE}, replace 0 intensities with 1 before applying mmgmos.}
  \item{gsnorm}{character. specifying the algorithm of global scaling normalisation.}
  \item{savepar}{Logical value. If \code{TRUE} the estimated parameters of the model are saved in file par\_mmgmos.txt and phi\_mmgmos.txt. }
  \item{eps}{Optimisation termination criteria.}
  \item{orig.phis}{Logical value. If \code{TRUE}, use phi values created from hgu95a array.}
}
\details{
  The obtained expression measures are in log base 2 scale.
   
  The algorithms of global scaling normalisation can be one of "mean", "median" and "meanlog".
  "mean" and "meanlog" are mean-centered normalisation on raw scale and log scale respectively, and "median" 
  is median-centered normalisation.
 
  There are 2*n+2 columns in file par\_mmgmos.txt, n is the number of chips. The first n columns are 'alpha' values for n chips, the next n columns are 'a' values for n 
  chips, column 2*n+1 is 'c' values and the final column is values for 'd'. The file phi\_mmgmos.txt keeps the final optimal value of 'phi'.
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

\author{Xuejun Liu, Magnus Rattray, Marta Milo, Neil D. Lawrence}

\seealso{ Related class \code{\link{exprReslt-class}} and related method \code{\link[mmgmos]{mgmos}}}
\examples{
## load example data from package Affy
data(affybatch.example)

## use method mmgmos to calculate the expression levels and related confidence 
## of the measures for the example data
eset<-mmgmos(affybatch.example)
}
\keyword{manip}
