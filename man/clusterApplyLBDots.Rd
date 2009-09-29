\name{clusterApplyLBDots}
\alias{clusterApplyLBDots}
\title{clusterApplyLB with dots to indicate progress}
\description{
This is basically the \code{\link[snow:snow-cluster]{clusterApplyLB}}
function from the \pkg{snow} package, but with dots displayed to
indicate progress. 
}
\usage{
clusterApplyLBDots(cl, x, fun, ...)
}
\arguments{
	\item{cl}{ cluster object }
	\item{x}{ array }
	\item{fun}{ function or character string naming a function }
	\item{\dots}{ additional arguments to pass to standard function }
}
\author{ Richard D. Pearson (modified from original \pkg{snow} function)}
\keyword{manip}
