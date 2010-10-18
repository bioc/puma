\name{pumaPCAExpectations-class}
\docType{class}
\alias{pumaPCAExpectations-class}
\alias{pumaPCAExpectations}

\title{Class pumaPCAExpectations}
\description{
This is a class representation for storing a set of expectations from a pumaPCA model. It is an internal representation and shouldn't normally be instantiated.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("pumaPCAExpectations", ...)}.
}
\section{Slots}{
	\describe{
		\item{\code{x}:}{Object of class "matrix" representing x}
		\item{\code{xxT}:}{Object of class "array" representing xxT}
		\item{\code{logDetCov}:}{Object of class "numeric" representing logDetCov}
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
