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
A separate list will be created for each contrast of interest.

While it is possible to run this function on data from individual arrays, it is generally recommended that this function is run on the output of the function \code{\link{pumaComb}} (which combines information from replicates).
}
\value{
	An object of class \code{\link{DEResult-class}}.
}
\author{ Richard D. Pearson }
\seealso{Related methods \code{\link{calculateLimma}}, \code{\link{calculateFC}}, \code{\link{calculateTtest}}, \code{\link{pumaComb}}, \code{\link{mmgmos}}, \code{\link{pplr}}, \code{\link{createDesignMatrix}} and \code{\link{createContrastMatrix}}}
\examples{
	library(puma)
	data(affybatch.example)
	pData(affybatch.example) <- data.frame("level"=c("twenty","twenty","ten")
	    , "batch"=c("A","B","A"), row.names=rownames(pData(affybatch.example)))
	eset_mmgmos <- mmgmos(affybatch.example)
	eset_mmgmos_normd <- pumaNormalize(eset_mmgmos)
	eset_comb <- pumaComb(eset_mmgmos_normd)
	esetDE <- pumaDE(eset_comb)
	esetDE$genes[1:6,]
	esetDE$p[1:6,]
}
\keyword{manip}
