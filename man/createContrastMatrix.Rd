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
	\item{eset}{ An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. }
	\item{design}{ A design matrix }
%	\item{others}{ A Boolean specifying whether to create vs\_other contrasts. }
}
\details{
The \pkg{puma} package has been designed to be as easy to use as possible, while not compromising on power and flexibility. One of the most difficult tasks for many users, particularly those new to microarray analysis, or statistical analysis in general, is setting up design and contrast matrices. The \pkg{puma} package will automatically create such matrices, and we believe the way this is done will suffice for most users' needs.

It is important to recognise that the automatic creation of design and contrast matrices will only happen if appropriate information about the levels of each factor is available for each array in the experimental design. This data should be held in an \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}} class. The easiest way of doing this is to ensure that the \code{\link[Biobase:class.AnnotatedDataFrame]{AnnotatedDataFrame}} object holding the raw CEL file data has an appropriate \code{phenoData} slot. This information will then be passed through to any \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object created, for example through the use of \code{\link{mmgmos}}. The \code{phenoData} slot of an \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}} object can also be manipulated directly if necessary.

Design and contrast matrices are dependent on the experimental design. The simplest experimental designs have just one factor, and hence the \code{phenoData} slot will have a matrix with just one column. In this case, each unique value in that column will be treated as a distinct level of the factor, and hence \code{\link{pumaComb}} will group arrays according to these levels. If there are just two levels of the factor, e.g. A and B, the contrast matrix will also be very simple, with the only contrast of interest being A vs B. For factors with more than two levels, a contrast matrix will be created which reflects all possible combinations of levels. For example, if we have three levels A, B and C, the contrasts of interest will be A vs B, A vs C and B vs C. In addition, if the \code{others} argument is set to TRUE, the following additional contrasts will be created: A vs other (i.e. A vs B \& C), B vs other and C vs other. Note that these additional contrasts are experimental, and not currently recommended for use in calculating differential expression.

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
# This is a simple example based on a real data set. Note that this is an "unbalanced" design, the "level" factor has two replicates of the "twenty" condition, but only one replicate of the "ten" condition. Also note that the second factor, "batch" is not used in the design or contrast matrices, as we don't have every combination of the levels of "level" and "batch" (there is no array for level=twenty and batch=B).

#	Next 4 lines commented out to save time in package checks, and saved version used
# if (require(affydata)) {
#	data(Dilution)
#	eset_mmgmos <- mmgmos(Dilution)
# }
data(eset_mmgmos)
createContrastMatrix(eset_mmgmos)

# The following shows a set of 15 synthetic data sets with increasing complexity. We first create the data sets, then look at the contrast matrices.

# single 2-level factor
eset1 <- new("ExpressionSet", exprs=matrix(0,100,4))
pData(eset1) <- data.frame("class"=c(1,1,2,2))

# single 2-level factor - unbalanced design
eset2 <- new("ExpressionSet", exprs=matrix(0,100,4))
pData(eset2) <- data.frame("class"=c(1,2,2,2))

# single 3-level factor
eset3 <- new("ExpressionSet", exprs=matrix(0,100,6))
pData(eset3) <- data.frame("class"=c(1,1,2,2,3,3))

# single 4-level factor
eset4 <- new("ExpressionSet", exprs=matrix(0,100,8))
pData(eset4) <- data.frame("class"=c(1,1,2,2,3,3,4,4))

# 2x2 factorial
eset5 <- new("ExpressionSet", exprs=matrix(0,100,8))
pData(eset5) <- data.frame("fac1"=c("a","a","a","a","b","b","b","b"), "fac2"=c(1,1,2,2,1,1,2,2))

# 2x2 factorial - unbalanced design
eset6 <- new("ExpressionSet", exprs=matrix(0,100,10))
pData(eset6) <- data.frame("fac1"=c("a","a","a","b","b","b","b","b","b","b"), "fac2"=c(1,2,2,1,1,1,2,2,2,2))

# 3x2 factorial
eset7 <- new("ExpressionSet", exprs=matrix(0,100,12))
pData(eset7) <- data.frame("fac1"=c("a","a","a","a","b","b","b","b","c","c","c","c"), "fac2"=c(1,1,2,2,1,1,2,2,1,1,2,2))

# 2x3 factorial
eset8 <- new("ExpressionSet", exprs=matrix(0,100,12))
pData(eset8) <- data.frame(
	"fac1"=c("a","a","a","a","a","a","b","b","b","b","b","b")
,	"fac2"=c(1,1,2,2,3,3,1,1,2,2,3,3) )

# 2x2x2 factorial
eset9 <- new("ExpressionSet", exprs=matrix(0,100,8))
pData(eset9) <- data.frame(
	"fac1"=c("a","a","a","a","b","b","b","b")
,	"fac2"=c(1,1,2,2,1,1,2,2)
,	"fac3"=c("X","Y","X","Y","X","Y","X","Y") )

# 3x2x2 factorial
eset10 <- new("ExpressionSet", exprs=matrix(0,100,12))
pData(eset10) <- data.frame(
	"fac1"=c("a","a","a","a","b","b","b","b","c","c","c","c")
,	"fac2"=c(1,1,2,2,1,1,2,2,1,1,2,2)
,	"fac3"=c("X","Y","X","Y","X","Y","X","Y","X","Y","X","Y") )

# 3x2x2 factorial
eset11 <- new("ExpressionSet", exprs=matrix(0,100,12))
pData(eset11) <- data.frame(
	"fac1"=c("a","a","a","a","a","a","b","b","b","b","b","b")
,	"fac2"=c(1,1,2,2,3,3,1,1,2,2,3,3)
,	"fac3"=c("X","Y","X","Y","X","Y","X","Y","X","Y","X","Y") )

# 3x2x2 factorial
eset12 <- new("ExpressionSet", exprs=matrix(0,100,18))
pData(eset12) <- data.frame(
	"fac1"=c("a","a","a","a","a","a","b","b","b","b","b","b","c","c","c","c","c","c")
,	"fac2"=c(1,1,2,2,3,3,1,1,2,2,3,3,1,1,2,2,3,3)
,	"fac3"=c("X","Y","X","Y","X","Y","X","Y","X","Y","X","Y","X","Y","X","Y","X","Y") )

# 2x2x2x2 factorial
eset13 <- new("ExpressionSet", exprs=matrix(0,100,16))
pData(eset13) <- data.frame(
	"fac1"=c("a","a","a","a","a","a","a","a","b","b","b","b","b","b","b","b")
,	"fac2"=c(0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1)
,	"fac3"=c(2,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3)
,	"fac4"=c("X","Y","X","Y","X","Y","X","Y","X","Y","X","Y","X","Y","X","Y") )

# "Un-analysable" data set - all arrays are from the same class
eset14 <- new("ExpressionSet", exprs=matrix(0,100,4))
pData(eset14) <- data.frame("class"=c(1,1,1,1))

# "Non-factorial" data set - there are no arrays for fac1="b" and fac2=2. In this case only the first factor (fac1) is used.
eset15 <- new("ExpressionSet", exprs=matrix(0,100,6))
pData(eset15) <- data.frame("fac1"=c("a","a","a","a","b","b"), "fac2"=c(1,1,2,2,1,1))

createContrastMatrix(eset1)
createContrastMatrix(eset2)
createContrastMatrix(eset3)
createContrastMatrix(eset4)
createContrastMatrix(eset5)
createContrastMatrix(eset6)
createContrastMatrix(eset7)
createContrastMatrix(eset8)
createContrastMatrix(eset9)
# For the last 4 data sets, the contrast matrices get pretty big, so we'll just show the names of each contrast
colnames(createContrastMatrix(eset10))
colnames(createContrastMatrix(eset11))
# Note that the number of contrasts can rapidly get very large for multi-factorial experiments!
colnames(createContrastMatrix(eset12))
# For this final data set, note that the puma package does not currently create interaction terms for data sets with 4 of more factors. E-mail the author if you would like to do this.
colnames(createContrastMatrix(eset13))

# "Un-analysable" data set - all arrays are from the same class - gives an error.  Note that we've commented this out so that we don't get errors which would make the package fail the Bioconductor checks!
# createContrastMatrix(eset14)
# "Non-factorial" data set - there are no arrays for fac1="b" and fac2=2. In this case only the first factor (fac1) is used.
createContrastMatrix(eset15)

}
\keyword{manip}
