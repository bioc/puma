\name{DEResult}
\docType{class}
\alias{class:DEResult}
\alias{DEResult}

\alias{DEResult-class}

\alias{statistic}
\alias{statistic<-}
\alias{FC}
\alias{FC<-}
\alias{statisticDescription}
\alias{statisticDescription<-}
\alias{DEMethod}
\alias{DEMethod<-}
\alias{pLikeValues}
\alias{topGenes}
\alias{topGeneIDs}
\alias{numberOfProbesets}
\alias{numberOfGenes}
\alias{numberOfContrasts}

\alias{statistic,DEResult-method}
\alias{statistic<-,DEResult,matrix-method}
\alias{FC,DEResult-method}
\alias{FC<-,DEResult,matrix-method}
\alias{statisticDescription,DEResult-method}
\alias{statisticDescription<-,DEResult,character-method}
\alias{DEMethod,DEResult-method}
\alias{DEMethod<-,DEResult,character-method}
\alias{pLikeValues,DEResult-method}
\alias{topGenes,DEResult-method}
\alias{topGeneIDs,DEResult-method}
\alias{numberOfProbesets,DEResult-method}
\alias{numberOfGenes,DEResult-method}
\alias{numberOfContrasts,DEResult-method}
\alias{show,DEResult-method}
\alias{write.reslts,DEResult-method}

\title{
  Class DEResult
}

\description{
  Class to contain and describe results of a differential expression (DE) analysis. The main components are \code{statistic} which hold the results of any statistic (e.g. p-values, PPLR values, etc.), and \code{FC} which hold the fold changes.
}
\section{Creating Objects}{

  \code{DEResult} objects will generally be created using one of the functions \code{\link{pumaDE}}, \code{\link{calculateLimma}}, \code{\link{calculateFC}} or \code{\link{calculateTtest}}.

  Objects can also be created from scratch:

  \code{new("DEResult")}

  \code{new("DEResult",
		statistic=matrix()
	,	FC=matrix()
	,	statisticDescription="unknown"
	,	DEMethod="unknown"
	)
  }
}
\section{Slots}{
   \describe{
     \item{\code{statistic}:}{Object of class "matrix" holding the statistics returned by the DE method.}
     \item{\code{FC}:}{Object of class "matrix" holding the fold changes returned by the DE method.}
     \item{\code{statisticDescription}:}{A text description of the contents of the \code{statistic} slot.}
     \item{\code{DEMethod}:}{A string indicating which DE method was used to create the object.}
  }
}
\section{Methods}{

  Class-specific methods.
  \describe{
    \item{\code{statistic(DEResult)}, \code{statistic(DEResult,matrix)<-}}{Access and
      set the \code{statistic} slot.}
    \item{\code{FC(DEResult)}, \code{FC(DEResult,matrix)<-}}{Access and
      set the \code{FC} slot.}
    \item{\code{statisticDescription(DEResult)}, \code{statisticDescription(DEResult,character)<-}}{Access and
      set the \code{statisticDescription} slot.}
    \item{\code{DEMethod(DEResult)}, \code{DEMethod(DEResult,character)<-}}{Access and
      set the \code{DEMethod} slot.}
    \item{\code{pLikeValues(object, contrast=1, direction="either")}}{Access the statistics of an object of class \code{DEResult}, converted to "p-like values".  If the object holds information on more than one contrast, only the values of the statistic for contrast number \code{contrast} are given. Direction can be "either" (meaning we want order genes by probability of being either up- or down-regulated), "up" (meaning we want to order genes by probability of being up-regulated), or "down" (meaning we want to order genes by probability of being down-regulated). "p-like values" are defined as values between 0 and 1, where 0 identifies the highest probability of being differentially expressed, and 1 identifies the lowest probability of being differentially expressed. We use this so that we can easily compare results from methods that provide true p-values (e.g. \code{\link{calculateLimma}}) and methods methods that do not provide p-values (e.g. \code{\link{pumaDE}}). For objects created using \code{\link{pumaDE}}, this returns 1-PPLR if the direction is "up", PPLR if direction is "down", and 1-abs(2*(PPLR-0.5)) if direction is "either". For objects created using \code{\link{calculateLimma}} or \code{\link{calculateTtest}}, this returns the p-value if direction is "either", ((p-1 * sign(FC))/2)+ 0.5, if the direction is "up", and ((1-p * sign(FC))/2)+ 0.5 if the direction is "down". For all other methods, this returns the rank of the appropriate statistic, scaled to lie between 0 and 1. \code{contrast} will be returned.}
    \item{\code{topGenes(object, numberOfGenes=1, contrast=1, direction="either")}}{Returns the index numbers (row numbers) of the genes determined to be most likely to be differentially expressed. \code{numberOfGenes} specifies the number of genes to be returned by the function. If the object holds information on more than one contrast, only the values of the statistic for contrast number \code{contrast} are given. Direction can be "either" (meaning we want order genes by probability of being either up- or down-regulated), "up" (meaning we want to order genes by probability of being up-ragulated), or "down" (meaning we want to order genes by probability of being down-regulated). Note that genes are ordered by "p-like values" (see \code{pLikeValues}). \code{object} is an object of class \code{DEResult}.}
    \item{\code{topGeneIDs(object, numberOfGenes=1, contrast=1, direction="either")}}{Returns the Affy IDs (row names) of the genes determined to be most likely to be differentially expressed. \code{numberOfGenes} specifies the number of genes to be returned by the function. If the object holds information on more than one contrast, only the values of the statistic for contrast number \code{contrast} are given. Direction can be "either" (meaning we want order genes by probability of being either up- or down-regulated), "up" (meaning we want to order genes by probability of being up-ragulated), or "down" (meaning we want to order genes by probability of being down-regulated). Note that genes are ordered by "p-like values" (see \code{pLikeValues}). \code{object} is an object of class \code{DEResult}.}
    \item{\code{numberOfProbesets(object)}}{Returns the number of probesets (number of rows) in an object of class \code{DEResult}. This method is synonymous with numberOfGenes.}
    \item{\code{numberOfGenes(object)}}{Returns the number of probesets (number of rows) in an object of class \code{DEResult}. This method is synonymous with numberOfProbesets.}
    \item{\code{numberOfContrasts(object)}}{Returns the number of contrasts (number of columns) in an object of class \code{DEResult}.}
    \item{\code{write.reslts(object)}}{\code{signature(x = "DEResult")}: writes the statistics and related fold changes (FCs) to
      files. It takes the same arguments as \code{\link{write.table}}. The argument "file" does not need to set any 
      extension. The different file marks and extension "csv" will be added automatically. The default file name is "tmp". 
      In the final results, statistics are in the file "tmp\_statistics.csv", and FCs are in 
      "tmp\_FCs.csv" respectively. }	}

  Standard generic methods:
	\describe{

		\item{\code{show(object)}}{Informatively display object contents.}
	}
}

\author{ Richard D. Pearson }

\seealso{
  Related methods \code{\link{pumaDE}}, \code{\link{calculateLimma}}, \code{\link{calculateFC}} or \code{\link{calculateTtest}}.
}
\examples{
## Create an example DEResult object
#	Next 4 lines commented out to save time in package checks, and saved version used
# if (require(affydata)) {
#	data(Dilution)
#	eset_mmgmos <- mmgmos(Dilution)
# }
data(eset_mmgmos)

#	Next line used so eset_mmgmos only has information about the liver factor
#	The scanner factor will thus be ignored, and the two arrays of each level
#	of the liver factor will be treated as replicates
pData(eset_mmgmos) <- pData(eset_mmgmos)[,1,drop=FALSE]

#	To save time we'll just use 100 probe sets for the example
eset_mmgmos_100 <- eset_mmgmos[1:100,]
eset_comb <- pumaComb(eset_mmgmos_100)

esetDE <- pumaDE(eset_comb)

## Use some of the methods
statisticDescription(esetDE)
DEMethod(esetDE)
numberOfProbesets(esetDE)
numberOfContrasts(esetDE)
topGenes(esetDE)
topGenes(esetDE, 3)
pLikeValues(esetDE)[topGenes(esetDE,3)]
topGeneIDs(esetDE, 3)
topGeneIDs(esetDE, 3, direction="down")

## save the expression results into files
write.reslts(esetDE, file="example")
}
\keyword{classes}
