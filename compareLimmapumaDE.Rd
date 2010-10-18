\name{compareLimmapumaDE}
\alias{compareLimmapumaDE}
\title{Compare pumaDE with a default Limma model}
\description{
This function compares the identification of differentially expressed (DE) genes using the \code{\link{pumaDE}} function and the \pkg{limma} package.
}
\usage{
compareLimmapumaDE(
	eset_mmgmos
,	eset_comb = NULL
,	eset_other = eset_mmgmos
,	limmaRes = calculateLimma(eset_other)
,	pumaDERes = pumaDE(eset_comb)
,	contrastMatrix = createContrastMatrix(eset_mmgmos)
,	numberToCompareForContrasts = 3
,	numberToCompareForVenn = 100
,	plotContrasts = TRUE
,	contrastsFilename = NULL
,	plotOther = FALSE
,	otherFilename = "other"
,	plotBcombContrasts = FALSE
,	bcombContrastsFilename = "bcomb_contrasts"
,	plotVenn = FALSE
,	vennFilename = "venn.pdf"
,	showTopMatches = FALSE
,	returnResults = FALSE
)
}
\arguments{
	\item{eset_mmgmos}{ An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}, that includes both expression levels as well as standard errors of the expression levels. This will often have been created using \code{\link{mmgmos}}, but might also have been created by \code{\link{mgmos}}, or any other method capable of providing standard errors. }
	\item{eset_comb}{ An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}, includes both expression levels as well as standard errors of the expression levels for each unique condition in an experiment (i.e. created from combining the information from each replicate). This will usually have been created using  \code{\link{pumaComb}}. }
	\item{eset_other}{ An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}, that includes expression levels , and may optionally also include standard errors of the expression levels. This is used for comparison with \code{eset_mmgmos}, and might have been created by any summarisation method, e.g. \code{\link[affy]{rma}}. }
	\item{limmaRes}{ A list with two elements, usually created using the function \code{\link{calculateLimma}}. The first element is a matrix of p-values. Each column represent one contrast. Within each column the p-values are ordered. The second element is a matrix of row numbers, which can be used to map p-values back to probe sets. If not supplied this will be automatically created from \code{eset_other}. }
	\item{pumaDERes}{ A list with two elements, usually created using the function \code{\link{pumaDE}}. The first element is a matrix of PPLR values. Each column represent one contrast. Within each column the PPLR values are ordered. The second element is a matrix of row numbers, which can be used to map PPLR values back to probe sets. If not supplied this will be automatically created from \code{eset_comb}. }
	\item{contrastMatrix}{ A contrast matrix. If not supplied this will be created from \code{eset_mmgmos} }
	\item{numberToCompareForContrasts}{ An integer specifying the number of most differentially expressed probe sets (genes) that will be used in comparison charts. }
	\item{numberToCompareForVenn}{ An integer specifying the number of most differentially expressed probe sets (genes) that will be used for comparison in the Venn  diagram. }
	\item{plotContrasts}{ A boolean specifying whether or not to plot the most differentially expressed probe sets (genes) for each contrast for the \code{eset_mmgmos} \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. }
	\item{contrastsFilename}{ A character string specifying a file name stem for the PDF files which will be created to hold the contrast plots for the \code{eset_mmgmos} \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. The actually filenames will have the name of the contrast appended to this stem. }
	\item{plotOther}{ A boolean specifying whether or not to plot the most differentially expressed probe sets (genes) for each contrast for the \code{eset_other} \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. }
	\item{otherFilename}{ A character string specifying a file name stem for the PDF files which will be created to hold the contrast plots for the \code{eset_other} \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. The actually filenames will have the name of the contrast appended to this stem. }
	\item{plotBcombContrasts}{ A boolean specifying whether or not to plot the most differentially expressed probe sets (genes) for each contrast for the \code{eset_comb} \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. }
	\item{bcombContrastsFilename}{ A character string specifying a file name stem for the PDF files which will be created to hold the contrast plots for the \code{eset_comb} \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}. The actually filenames will have the name of the contrast appended to this stem. }
	\item{plotVenn}{ A boolean specifying whether or not to plot a Venn diagram showing the overlap in the most differentially expressed probe sets (genes) as identified from the two different methods being compared. }
	\item{vennFilename}{ A character string specifying the filename for the PDF file which will hold the Venn diagram showing the overlap in the most differentially expressed probe sets (genes) as identified from the two different methods being compared. }
	\item{showTopMatches}{ A boolean specifying whether or not to show the probe sets which are deemed most likely to be differentially expressed. }
	\item{returnResults}{ A boolean specifying whether or not to return a list containing results generated. }
}
\value{
The main outputs from this function are a number of PDF files.

The function only returns results if returnResults=TRUE
}
\author{ Richard D. Pearson }
\seealso{Related methods \code{\link{pumaDE}} and \code{\link{calculateLimma}}}
\keyword{manip}
