\name{create_eset_r}
\alias{create_eset_r}
\title{Create an ExpressionSet from a PPLR matrix}
\description{
This is really an internal function called from pumaComb. It is used to create an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object from the output of the bcomb function (which was originally part of the \pkg{pplr} package. Don't worry about it!
}
\usage{
create_eset_r(
	eset
,	r
,	design.matrix=createDesignMatrix(eset)
)
}
\arguments{
	\item{eset}{ An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. The phenotype information from this is used as the phenotype information of the returned object }
	\item{r}{ A data frame with components named 'M1', 'M2', and so on, which represent the mean expression values for condition 1, condition 2, and so on. It also has components named 'Std1', 'Std2', and so on, which represent the standard deviation of the gene expression values for condition 1, condtion 2, and so on. This type of data frame is output by function \code{\link{bcomb}} and \code{\link{hcomb}} }
	\item{design.matrix}{ A design matrix. }
}
\value{
An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.
}
\author{ Richard D. Pearson }
\seealso{Related methods \code{\link{bcomb}}, \code{\link{hcomb}}, \code{\link{pumaComb}} and \code{\link{pumaCombImproved}}}
\keyword{manip}
