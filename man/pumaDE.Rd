\name{pumaDE}
\alias{pumaDE}
\title{Calculate differential expression between conditions}
\description{
The function generates lists of genes ranked by probability of differential expression (DE). This uses the PPLR method.
}
\usage{
pumaDE(
	eset
,	design.matrix = createDesignMatrix(eset)
,	contrast.matrix = createContrastMatrix(eset)
)
}
\arguments{
	\item{eset}{ An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. }
	\item{design.matrix}{ A design matrix }
	\item{contrast.matrix}{ A contrast matrix }
}
\details{
A separate list of genes will be created for each contrast of interest.

Note that this class returns a \code{\link{DEResult-class}} object. This object contains information on both the PPLR statistic values (which should generally be used to rank genes for differential expression), as well as fold change values (which are generally not recommended for ranking genes, but which might be useful, for example, to use as a filter). To understand more about the object returned see \code{\link{DEResult-class}}, noting that when created a DEResult object with the pumaDE function, the \code{statistic} method should be used to return PPLR values. Also note that the \code{pLikeValues} method can be used on the returned object to create values which can more readily be compared with p-values returned by other methods such as variants of t-tests (limma, etc.).

While it is possible to run this function on data from individual arrays, it is generally recommended that this function is run on the output of the function \code{\link{pumaComb}} (which combines information from replicates).
}
\value{
	An object of class \code{\link{DEResult-class}}.
}
\author{ Richard D. Pearson }
\seealso{Related methods \code{\link{calculateLimma}}, \code{\link{calculateFC}}, \code{\link{calculateTtest}}, \code{\link{pumaComb}}, \code{\link{mmgmos}}, \code{\link{pplr}}, \code{\link{createDesignMatrix}} and \code{\link{createContrastMatrix}}}
\examples{
	data(affybatch.example)
	varLabels(affybatch.example) <- c("level")
	pData(affybatch.example) <- data.frame("level"=c("twenty","twenty","ten")
	    , "batch"=c("A","B","A"), row.names=rownames(pData(affybatch.example)))
	eset_mmgmos <- mmgmos(affybatch.example)
	eset_comb <- pumaComb(eset_mmgmos)
	pumaDEResults <- pumaDE(eset_comb)

	topGeneIDs(pumaDEResults,6) # Gives probeset identifiers
	topGenes(pumaDEResults,6) # Gives row numbers
	statistic(pumaDEResults)[topGenes(pumaDEResults,6),] # PPLR scores of top six genes
	FC(pumaDEResults)[topGenes(pumaDEResults,6),] # Fold-change of top six genes
}
\keyword{manip}
