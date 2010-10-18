\name{pumaPCARes-class}
\docType{class}
\alias{class:pumaPCARes}
\alias{pumaPCARes}
\alias{pumaPCARes-class}
\alias{write.reslts,pumaPCARes-method}

\title{Class pumaPCARes}
\description{
This is a class representation for storing the outputs of the pumaPCA function. Objects of this class should usually only be created through the \code{\link{pumaPCA}} function.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("pumaPCARes", ...)}.
}
\section{Slots}{
	\describe{
		\item{\code{model}:}{Object of class "pumaPCAModel" representing the model parameters}
		\item{\code{expectations}:}{Object of class "pumaPCAExpectations" representing the model expectations}
		\item{\code{varY}:}{Object of class "matrix" representing the variance in the expression levels}
		\item{\code{Y}:}{Object of class "matrix" representing the expression levels}
		\item{\code{phenoData}:}{Object of class "AnnotatedDataFrame" representing the phenotype information}
		\item{\code{timeToCompute}:}{Object of class "numeric" representing the time it took \code{\link{pumaPCA}} to run}
		\item{\code{numberOfIterations}:}{Object of class "numeric" representing the number of iterations it took \code{\link{pumaPCA}} to converge}
		\item{\code{likelihoodHistory}:}{Object of class "list" representing the history of likelihood values while \code{\link{pumaPCA}} was running}
		\item{\code{timingHistory}:}{Object of class "list" representing the history of how long each iteration took while \code{\link{pumaPCA}} was running}
		\item{\code{modelHistory}:}{Object of class "list" representing the history of how the model was changing while \code{\link{pumaPCA}} was running}
		\item{\code{exitReason}:}{Object of class "character" representing the reason \code{\link{pumaPCA}} halted. Can take the values "Update of Likelihood less than tolerance x", "Update of W less than tolerance x", "Iterations exceeded", "User interrupt", "unknown exit reason"}
	}
}
\section{Methods}{
	\describe{
		\item{plot}{\code{signature(x="pumaPCARes-class")}: plots two principal components on a scatter plot.}
		\item{write.reslts}{\code{signature(x = "pumaPCARes-class")}: writes the principal components for each array to a file. It takes the same arguments as \code{\link{write.table}}. The argument "file" does not need to set any extension. The file name and extension "csv" will be added automatically. The default file name is "tmp". }
	}
}
\author{ Richard D. Pearson }

\seealso{
Related method \code{\link{pumaPCA}} and related class \code{\link{pumaPCARes}}.
}
\keyword{classes}
