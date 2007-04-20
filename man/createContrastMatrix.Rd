\name{createContrastMatrix}
\alias{createContrastMatrix}
\title{Automatically create a contrast matrix from an ExpressionSet and optional design matrix}
\description{
To appear
}
\usage{
createContrastMatrix(eset, design=NULL)
}
\arguments{
	\item{eset}{ An object of class \code{\link[Biobase]{ExpressionSet}}. }
	\item{design}{ A design matrix }
}
\details{
The \pkg{puma} package has been designed to be as easy to use as possible, while not compromising on power and flexibility. One of the most difficult tasks for many users, particularly those new to microarray analysis, or statistical analysis in general, is setting up design and contrast matrices. The \pkg{puma} package will automatically create such matrices, and we believe the way this is done will suffice for most users' needs.

It is important to recognise that the automatic creation of design and contrast matrices will only happen if appropriate information about the levels of each factor is available for each array in the experimental design. This data should be held in an \code{\link[Biobase]{AnnotatedDataFrame}} class. The easiest way of doing this is to ensure that the \code{\link[Biobase]{AnnotatedDataFrame}} object holding the raw CEL file data has an appropriate \code{phenoData} slot. This information will then be passed through to any \code{\link[Biobase]{ExpressionSet}} object created, for example through the use of \code{\link{mmgmos}}. The \code{phenoData} slot of an \code{\link[Biobase]{ExpressionSet}} object can also be manipulated directly if necessary.

Design and contrast matrices are dependent on the experimental design. The simplest experimental designs have just one factor, and hence the \code{phenoData} slot will have a matrix with just one column. In this case, each unique value in that column will be treated as a distinct level of the factor, and hence \code{\link{pumaComb}} will group arrays according to these levels. If there are just two levels of the factor, e.g. A and B, the contrast matrix will also be very simple, with the only contrast of interest being A vs B. For factors with more than two levels, a contrast matrix will be created which reflects all possible combinations of levels. For example, if we have three levels A, B and C, the contrasts of interest will be A vs B, A vs C and B vs C.

If we now consider the case of two or more factors, things become more complicated. There are now two cases to be considered: factorial experiments, and non-factorial experiments. A factorial experiment is one where all the combinations of the levels of each factor are tested by at least one array (though ideally we would have a number of biological replicates for each combination of factor levels). The \code{estrogen} case study from the package vignette is an example of a factorial experiment.

A non-factorial experiment is one where at least one combination of levels is not tested. If we treat the example used in the \code{\link{puma-package}} help page as a two-factor experiment (with factors ``level'' and ``batch''), we can see that this is not a factorial experiment as we have no array to test the conditions ``level=ten'' and ``batch=B''. We will treat the factorial and non-factorial cases separately in the following sections.

Factorial experiments

For factorial experiments, the design matrix will use all columns from the \code{phenoData} slot. This will mean that \code{\link{pumaComb}} will group arrays according to a combination of the levels of all the factors.

Non-factorial designs

For non-factorial designed experiments, we will simply ignore columns (right to left) from the \code{phenoData} slot until we have a factorial design or a single factor. We can see this in the example used in the \code{\link{puma-package}} help page. Here we have ignored the ``batch'' factor, and modelled the experiment as a single-factor experiment (with that single factor being ``level'').

}
\value{
  The result is a matrix. See the code below for an example.
}
\author{ Richard D. Pearson }
\seealso{
	Related methods \code{\link{createDesignMatrix}} and \code{\link{pumaDE}}
}
\examples{
	data(affybatch.example)
	pData(affybatch.example) <- data.frame("level"=c("twenty","twenty","ten")
	    , "batch"=c("A","B","A"), row.names=rownames(pData(affybatch.example)))
	eset_mmgmos <- mmgmos(affybatch.example)
	createContrastMatrix(eset_mmgmos)
}
\keyword{manip}
