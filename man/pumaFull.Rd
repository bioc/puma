\name{pumaFull}
\alias{pumaFull}
\title{Perform a full PUMA analysis}
\description{
Full analysis including pumaPCA and mmgmos/pumaDE vs rma/limma comparison
}
\usage{
pumaFull (
	affybatch = NULL
,	data_dir = getwd()
,	load_affybatch = FALSE
,	calculate_eset = TRUE
,	calculate_pumaPCAs = TRUE
,	calculate_bcomb = TRUE
,	mmgmosComparisons = FALSE
)

}
\arguments{
	\item{affybatch}{ ~~Describe \code{affybatch} here~~ }
	\item{data_dir}{ ~~Describe \code{data_dir} here~~ }
	\item{load_affybatch}{ ~~Describe \code{load_affybatch} here~~ }
	\item{calculate_eset}{ ~~Describe \code{calculate_eset} here~~ }
	\item{calculate_pumaPCAs}{ ~~Describe \code{calculate_pumaPCAs} here~~ }
	\item{calculate_bcomb}{ ~~Describe \code{calculate_bcomb} here~~ }
	\item{mmgmosComparisons}{ ~~Describe \code{mmgmosComparisons} here~~ }
}
\value{
  No return values. Various objects are saved as .rda files during the execution of this function, and various PDF files are created.
}
\author{ Richard D. Pearson }
\seealso{Related methods \code{\link{pumaDE}}, \code{\link{createDesignMatrix}} and \code{\link{createContrastMatrix}}}
\examples{
	data(affybatch.example)
	pData(affybatch.example) <- data.frame("level"=c("twenty","twenty","ten")
	    , "batch"=c("A","B","A"), row.names=rownames(pData(affybatch.example)))
	pumaFull(affybatch.example)
}
\keyword{manip}
