\name{gmhta}
\alias{gmhta}
\title{ Compute gene and transcript expression values and standard deviatons from hta2.0 CEL Files }
\description{
   This function converts an object of \code{\link[oligoClasses:FeatureSet-class]{FeatureSet}} into an object of class 
  \code{\link{exprReslt}} using the gamma model for hta2.0 chips. 
   This function obtains confidence of measures, standard deviation and 5, 
  25, 50, 75 and 95 percentiles, as well as the estimated expression levels.
}
\usage{
gmhta(
      object
     ,background=FALSE
     ,gsnorm=c("median", "none", "mean", "meanlog")
     ,savepar=FALSE
     ,eps=1.0e-6
     ,addConstant = 0
     ,cl=NULL
     ,BatchFold=10
)	
}
\arguments{
	\item{object}{an object of \code{\link[oligoClasses:FeatureSet-class]{FeatureSet}}}
          \item{background}{Logical value. If \code{TRUE}, perform background correction before applying gmhta.}
          \item{gsnorm}{character. specifying the algorithm of global scaling normalisation.}
	\item{savepar}{Logical value. If \code{TRUE} the estimated parameters of the model are saved in file par\_gmhta.txt }
	\item{eps}{Optimisation termination criteria.}
          \item{addConstant}{numeric. This is an experimental feature and should not generally be changed from the default value. }
          \item{cl}{This function can be parallelised by setting parameter cl. For more details, please refer to the vignette.}
          \item{BatchFold}{ we divide tasks into BatchFold*n jobs where n is the number of cluster nodes. The first n jobs are placed on the n nodes. When the first job is completed,the next job is placed on the available node. This continues until all jobs are completed. The default value is ten. The user also can change the value according to the number of cluster nodes n. We suggest that for bigger n BatchFold should be smaller.}
}
\details{
  The obtained expression measures are in log base 2 scale.
  Using the known relationships between genes, transcripts and probes,
  we propose a gamma model for hta2.0  data to calculate transcript and gene expression levels.
  The algorithms of global scaling normalisation can be one of "median", "none", "mean", "meanlog".
  "mean" and "meanlog" are mean-centered normalisation on raw scale and log scale respectively, and "median" 
  is median-centered normalisation. "none" will result in no global scaling normalisation being applied.
  This function can be parallelised by setting parameter cl. For more details, please refer to the vignette.
 

  }
\value{
  A list of two object of class \code{exprReslt}.
}
\references{XueJun Liu. (2013) puma 3.0: improved uncertainty propagation methods for gene and transcript expression
analysis, BMC Bioinformatics, 14:39.


Manhong Dai, Pinglang Wang,Andrew D. Boyd. (2005) Evolving gene/transcript definitions significantly
alter the interpretation of GeneChip data,Nucleic Acid Research 33(20):e175.
}

\author{Xuejun Liu,WuJun Zhang,Zhenzhu gao, Magnus Rattray}

\seealso{ Related class \code{\link{exprReslt-class}} }

\examples{
## The following scripts show the use of the method.
#library(puma)
## load CEL files
# object<-read.celfiles("celnames")
#eset<-gmhta(object,gsnorm="none",cl=cl)

}

\keyword{manip}
