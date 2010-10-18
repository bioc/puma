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
\seealso{Related methods \code{\link{calculateLimma}}, \code{\link{calculateFC}}, \code{\link{calculateTtest}}, \code{\link{pumaComb}}, \code{\link{pumaCombImproved}}, \code{\link{mmgmos}}, \code{\link{pplr}}, \code{\link{createDesignMatrix}} and \code{\link{createContrastMatrix}}}
\examples{
	#	Next 4 lines commented out to save time in package checks, and saved version used
    # if (require(affydata)) {
	#	data(Dilution)
	#	eset_mmgmos <- mmgmos(Dilution)
	# }
	data(eset_mmgmos)

	#	Next line shows that eset_mmgmos has 4 arrays, each of which is a different
	#   condition (the experimental design is a 2x2 factorial, with both liver and
	#	scanner factors)
	pData(eset_mmgmos)
	
	#	Next line shows expression levels of first 3 probe sets
	exprs(eset_mmgmos)[1:3,]

	#	Next line used so eset_mmgmos only has information about the liver factor
	#	The scanner factor will thus be ignored, and the two arrays of each level
	#	of the liver factor will be treated as replicates
	pData(eset_mmgmos) <- pData(eset_mmgmos)[,1,drop=FALSE]

	#	To save time we'll just use 100 probe sets for the example
	eset_mmgmos_100 <- eset_mmgmos[1:100,]
	eset_comb <- pumaComb(eset_mmgmos_100)
        eset_combimproved <- pumaCombImproved(eset_mmgmos_100)

	pumaDEResults <- pumaDE(eset_comb)
        pumaDEResults_improved <- pumaDE(eset_combimproved)

	topGeneIDs(pumaDEResults,6) # Gives probeset identifiers
	topGeneIDs(pumaDEResults_improved,6)
	topGenes(pumaDEResults,6) # Gives row numbers
        topGenes(pumaDEResults_improved,6)
	statistic(pumaDEResults)[topGenes(pumaDEResults,6),] # PPLR scores of top six genes
	statistic(pumaDEResults_improved)[topGenes(pumaDEResults_improved,6),] 
	FC(pumaDEResults)[topGenes(pumaDEResults,6),] # Fold-change of top six genes
	FC(pumaDEResults_improved)[topGenes(pumaDEResults_improved,6),] 
}
\keyword{manip}
