\name{justmgMOS}
\alias{justmgMOS}
\alias{just.mgmos}
\title{Compute mgmos Directly from CEL Files}
\description{
  This function converts CEL files into an \code{\link{exprReslt}} using mgmos.
}
\usage{
justmgMOS(\dots, filenames=character(0),
         widget=getOption("BioC")$affy$use.widgets,
         compress=getOption("BioC")$affy$compress.cel,
         celfile.path=getwd(),
         sampleNames=NULL,
         phenoData=NULL,
         description=NULL,
         notes="",
         background=TRUE, gsnorm=c("median", "none", "mean", "meanlog"), savepar=FALSE, eps=1.0e-6)

just.mgmos(\dots, filenames=character(0),
          phenoData=new("AnnotatedDataFrame"),
          description=NULL,
          notes="",
          compress=getOption("BioC")$affy$compress.cel,
          background=TRUE, gsnorm=c("median", "none", "mean", "meanlog"), savepar=FALSE, eps=1.0e-6)
}
\arguments{
  \item{\dots}{file names separated by comma.}
  \item{filenames}{file names in a character vector.}
  \item{widget}{a logical specifying if widgets should be used.}
  \item{compress}{are the CEL files compressed?}
  \item{celfile.path}{a character denoting the path where cel files locate.}        
  \item{sampleNames}{a character vector of sample names to be used in
    the \code{\link[oligoClasses:FeatureSet-class]{FeatureSet}}.}
 \item{phenoData}{an
          \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}}
          object.} 
  \item{description}{a \code{\link[Biobase:class.MIAME]{MIAME}} object.}
  \item{notes}{notes.}
  \item{background}{Logical value. If \code{TRUE}, then perform background correction before applying mgmos.}
  \item{gsnorm}{character. specifying the algorithm of global scaling normalisation.}
  \item{savepar}{Logical value. If \code{TRUE}, then the estimated parameters of the model are saved in file par\_mgmos.txt and phi\_mgmos.txt.}
  \item{eps}{Optimisation termination criteria.}
}
\details{
  This method should require much less RAM than the conventional
  method of first creating an \code{\link[oligoClasses:FeatureSet-class]{FeatureSet}} and then running
  \code{\link{mgmos}}.
  
  Note that this expression measure is given to you in log base 2 scale. This differs from
  most of the other expression measure methods.

  The algorithms of global scaling normalisation can be one of "median", "none", "mean", "meanlog".
  "mean" and "meanlog" are mean-centered normalisation on raw scale and log scale respectively, and "median" 
  is median-centered normalisation. "none" will result in no global scaling normalisation being applied.
}
\value{
 An \code{exprReslt}.
}
\seealso{Related class \code{\link{exprReslt-class}} and related method \code{\link{mgmos}}}

\keyword{manip}
