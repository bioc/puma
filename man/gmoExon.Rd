\name{gmoExon}
\alias{gmoExon}
\title{ Compute gene and transcript expression values and standard deviatons from exon CEL Files }
\description{
   This function converts an object of class \code{\link[affy:AffyBatch-class]{AffyBatch}} into an object of class 
  \code{\link{exprReslt}} using the gamma model for exon chips. 
   This function obtains confidence of measures, standard deviation and 5, 
  25, 50, 75 and 95 percentiles, as well as the estimated expression levels.
}
\usage{
gmoExon(
      object
     ,exontype = c("Human", "Mouse", "Rat")
     ,GT = c("gene", "transcript")
     ,background=FALSE
     ,gsnorm=c("median", "none", "mean", "meanlog")
     ,savepar=FALSE
     ,eps=1.0e-6
     ,addConstant = 0
)	
}
\arguments{
	\item{object}{an object of \code{\link[affy:AffyBatch-class]{AffyBatch}}}
	\item{exontype}{character. specifying the type of exon chip.}
          \item{GT}{character. specifying gene or transcript expression values and standard deviations to be obtained.}
          \item{background}{Logical value. If \code{TRUE}, perform background correction before applying gmoExon.}
          \item{gsnorm}{character. specifying the algorithm of global scaling normalisation.}
	\item{savepar}{Logical value. If \code{TRUE} the estimated parameters of the model are saved in file par\_gmoExon.txt }
	\item{eps}{Optimisation termination criteria.}
          \item{addConstant}{numeric. This is an experimental feature and should not generally be changed from the default value. }
}
\details{
  The obtained expression measures are in log base 2 scale.
  Using the known relationships between genes, transcripts and probes,
  we propose a gamma model for exon array data to calculate transcript and gene expression levels.
  The algorithms of global scaling normalisation can be one of "median", "none", "mean", "meanlog".
  "mean" and "meanlog" are mean-centered normalisation on raw scale and log scale respectively, and "median" 
  is median-centered normalisation. "none" will result in no global scaling normalisation being applied.
 

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
http://plato.la.asu.edu/donlp2.html

Risueno A, Fontanillo C, Dinger ME, De Las Rivas J. GATExplorer: genomic and transcriptomic explorer; mapping expression
probes to gene loci, transcripts, exons and ncRNAs. BMC Bioinformatics.2010.}

\author{Xuejun Liu, Zhenzhu gao, Magnus Rattray, Marta Milo, Neil D. Lawrence}

\seealso{ Related class \code{\link{exprReslt-class}} }

\examples{
## The following scripts show the use of the method.
## load CEL files
# datadir<-getwd()
# exonFilenames<-c("SR20070419HEX01.CEL", "SR20070419HEX02.CEL","SR20070419HEX06.CEL","SR20070419HEX07.CEL)
#affybatch.exon<-ReadAffy(filenames=exonFilenames, celfile.path=getwd());
# pData(affybatch.exon)<-data.frame("estrogen"=c("Tumor","Tumor","Normal","Normal"),"time.h"=c("10","10","48","48"),
#    row.names=rownames(pData(affybatch.exon)));
#affybatch.exon@cdfName<-"exon.pmcdf"
## use method gmoExon to calculate the expression levels and related confidence 
## of the measures for the example data
#eset<-gmoExon(affybatch.exon)
}

\keyword{manip}
