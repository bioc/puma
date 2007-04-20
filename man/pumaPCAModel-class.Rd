\name{pumaPCAModel-class}
\docType{class}
\alias{pumaPCAModel-class}
\alias{pumaPCAModel}

\title{Class pumaPCAModel}
\description{
This is a class representation for storing a pumaPCA model. It is an internal representation and shouldn't normally be instantiated.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("pumaPCAModel", ...)}.
}
\section{Slots}{
	\describe{
		\item{\code{sigma}:}{Object of class "numeric" representing sigma}
		\item{\code{m}:}{Object of class "matrix" representing m}
		\item{\code{Cinv}:}{Object of class "matrix" representing Cinv}
		\item{\code{W}:}{Object of class "matrix" representing W}
		\item{\code{mu}:}{Object of class "matrix" representing mu}
	}
}
\section{Methods}{
This class has no methods defined
}
\author{ Richard D. Pearson }

\seealso{
Related method \code{\link{pumaPCA}} and related class \code{\link{pumaPCARes}}.
}
\keyword{classes}
