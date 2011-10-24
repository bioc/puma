\name{exprReslt-class}
\docType{class}
\alias{class:exprReslt}
\alias{exprReslt}
\alias{exprReslt-class}
\alias{se.exprs}
\alias{se.exprs<-}
\alias{prcfifty}
\alias{prcfifty<-}
\alias{prcfive}
\alias{prcfive<-}
\alias{prcninfive}
\alias{prcninfive<-}
\alias{prcsevfive}
\alias{prcsevfive<-}
\alias{prctwfive}
\alias{prctwfive<-}
\alias{write.reslts}
\alias{show,exprReslt-method}
\alias{se.exprs,exprReslt-method}
\alias{se.exprs<-,exprReslt-method}
\alias{prcfifty,exprReslt-method}
\alias{prcfifty<-,exprReslt-method}
\alias{prcfive,exprReslt-method}
\alias{prcfive<-,exprReslt-method}
\alias{prcninfive,exprReslt-method}
\alias{prcninfive<-,exprReslt-method}
\alias{prcsevfive,exprReslt-method}
\alias{prcsevfive<-,exprReslt-method}
\alias{prctwfive,exprReslt-method}
\alias{prctwfive<-,exprReslt-method}
\alias{write.reslts,exprReslt-method}
\alias{write.reslts,ExpressionSet-method}

\title{Class exprReslt}
\description{This is a class representation for Affymetrix GeneChip probe level data. 
  The main component are the intensities, estimated expression levels and the confidence 
  of expression levels from multiple arrays 
  of the same \code{CDF} type. In extends \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("exprReslt", ...)}.
}
\section{Slots}{
  \describe{
  \item{\code{prcfive}:}{Object of class "matrix" representing the 5 percentile of the 
    observed expression levels. This is a matrix with columns representing patients or cases and 
    rows representing genes. }
  \item{\code{prctwfive}:}{Object of class "matrix" representing the 25 percentile of the 
    observed expression levels. This is a matrix with columns representing patients or cases and 
    rows representing genes.}
  \item{\code{prcfifty}:}{Object of class "matrix" representing the 50 percentile of the 
    observed expression levels. This is a matrix with columns representing patients or cases and 
    rows representing genes.}
  \item{\code{prcsevfive}:}{Object of class "matrix" representing the 75 percentile of the 
    observed expression levels. This is a matrix with columns representing patients or cases and 
    rows representing genes.}
  \item{\code{prcninfive}:}{Object of class "matrix" representing the 95 percentile of the 
    observed expression levels. This is a matrix with columns representing patients or cases and 
    rows representing genes.}
  \item{\code{phenoData}:}{Object of class "phenoData" inherited from
    \code{ExpressionSet}.}
  \item{\code{annotation}:}{A character string identifying the
    annotation that may be used for the \code{ExpressionSet} instance.}
  }
}
\section{Extends}{
Class \code{"ExpressionSet"}, directly.
}
\section{Methods}{
  \describe{
  \item{se.exprs}{\code{signature(object = "exprReslt")}: obtains the standard error of the estimated 
    expression levels. }
  \item{se.exprs<-}{\code{signature(object = "exprReslt")}: replaces the standard error of the estimated 
    expression levels. }
    \item{prcfifty}{\code{signature(object = "exprReslt")}: obtains the 50 percentile of the estimated 
      expression levels. }
    \item{prcfifty<-}{\code{signature(object = "exprReslt")}: replaces the 50 percentile of the estimated 
      expression levels. }
    \item{prcfive}{\code{signature(object = "exprReslt")}: obtains the 5 percentile of the estimated 
      expression levels.  }
    \item{prcfive<-}{\code{signature(object = "exprReslt")}: replaces the 5 percentile of the estimated 
      expression levels.  }
    \item{prcninfive}{\code{signature(object = "exprReslt")}: obtains the 95 percentile of the estimated 
      expression levels.  }
    \item{prcninfive<-}{\code{signature(object = "exprReslt")}: replaces the 95 percentile of the estimated 
      expression levels.  }
    \item{prcsevfive}{\code{signature(object = "exprReslt")}: obtains the 75 percentile of the estimated 
      expression levels.  }
    \item{prcsevfive<-}{\code{signature(object = "exprReslt")}: replaces the 75 percentile of the estimated 
      expression levels.  }
    \item{prctwfive}{\code{signature(object = "exprReslt")}: obtains the 25 percentile of the estimated 
      expression levels.  }
    \item{prctwfive<-}{\code{signature(object = "exprReslt")}: replaces the 25 percentile of the estimated 
      expression levels.  }
    \item{show}{\code{signature(object = "exprReslt")}: renders information about the exprReslt in a concise 
      way on stdout. }
    \item{write.reslts}{\code{signature(x = "exprReslt")}: writes the expression levels and related confidences to
      files. It takes the same arguments as \code{\link{write.table}}. The argument "file" does not need to set any 
      extension. The different file marks and extension "csv" will be added automatically. The default file name is "tmp". 
      In the final results, expression levels are in the file "tmp\_exprs.csv", standard deviations in 
      "tmp\_se.csv", 5 percentiles in "tmp\_prctile5.csv", likewise, 25, 50, 75 and 95 percentiles in "tmp\_prctile25.csv", 
      "tmp\_prctile50.csv", "tmp\_prctile75.csv" and "tmp\_prctile95.csv" respectively. }
  }
}
\author{ Xuejun Liu, Magnus Rattray, Marta Milo, Neil D. Lawrence, Richard D. Pearson }

\seealso{
  Related method \code{\link{mmgmos}} and related class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.
}
\examples{
## load example data from package affydata
#	Next 4 lines commented out to save time in package checks, and saved version used
# if (require(affydata)) {
#	data(Dilution)
#	eset_mmgmos <- mmgmos(Dilution)
# }
data(eset_mmgmos)

## save the expression results into files
write.reslts(eset_mmgmos, file="example")
}
\keyword{classes}
