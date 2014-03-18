\name{igmoExon}
\alias{igmoExon}
\title{ Separately Compute gene and transcript expression values and standard deviatons from exon CEL Files by the conditions.}
\description{
    The principle of this function is as same as the function gmoExon.This function separately calculates gene expression values by the conditions and 
    then combined every condition's results, and normalises them finally.                            
}
\usage{
igmoExon(
      cel.path
     ,SampleNameTable
     ,exontype = c("Human", "Mouse", "Rat")
     ,background=FALSE
     ,gsnorm=c("median", "none", "mean", "meanlog")     
     ,savepar=FALSE
     ,eps=1.0e-6
     ,addConstant = 0
     ,condition=c("Yes","No")
     ,cl=NULL
     ,BatchFold=10
)	
}
\arguments{
	\item{cel.path}{The directory where you put the CEL files.}
          \item{SampleNameTable}{It is a tab-separated  table with two columns,ordered by "Celnames","Condition"}
	\item{exontype}{character. specifying the type of exon chip.}
          \item{background}{Logical value. If \code{TRUE}, perform background correction before applying gmoExon.}
	\item{savepar}{Logical value. If \code{TRUE} the estimated parameters of the model are saved in file par\_gmoExon.txt }
	\item{eps}{Optimisation termination criteria.}
          \item{addConstant}{numeric. This is an experimental feature and should not generally be changed from the default value. }
          \item{gsnorm}{character. specifying the algorithm of global scaling normalisation.}
          \item{condition}{ Yes or No. ``Yes'' means  the \code{\link{igmoExon}} function separately calculates gene expression values by the conditions and then combined every condition's results, and normalises them finally. ``No'' means the \code{\link{igmoExon}} calulates the gene expression values as same as the \code{\link{gmoExon}} function.}
      \item{cl}{This function can be parallelised by setting parameter cl. For more details, please refer to the vignette.}
        \item{BatchFold}{ we divide tasks into BatchFold*n jobs where n is the number of cluster nodes. The first n jobs are placed on the n nodes. When the first job is completed,the next job is placed on the available node. This continues until all jobs are completed. The default value is ten. The user also can change the value according to the number of cluster nodes n. We suggest that for bigger n BatchFold should be smaller.}
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
  A list of two object of class \code{exprReslt}.
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
# cel.path<-cel.path;
# SampleNameTable<-"SampleNameTable"
#eset_igmoExon<-igmoExon(cel.path="cel.path"
                   # , SampleNameTable="SampleNameTable"
                   # , exontype="Human"
                   # , gsnorm="none", condition="Yes",cl=cl)
}
\keyword{manip}
