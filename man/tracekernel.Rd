\name{tracekernel}

\alias{tracekernel}
\title{Trace of product kernel smoother}


\description{
Evaluate the trace of the product of kernel smoother (Gaussian, Epanechnikov,
Uniform, Quartic). This function is not intended to be used directly.
}

\usage{
tracekernel(X,bx,kernelx,n,p)
}

\arguments{
  \item{X}{Matrix of explanatory variables, size \emph{n}, \emph{p}.}
  \item{bx}{The vector of different bandwidths, length \eqn{p}.}
  \item{kernelx}{Character string which allows to choose between gaussian kernel
  (\code{"g"}), Epanechnikov (\code{"e"}), uniform (\code{"u"}),
  quartic (\code{"q"}).}
  \item{n}{Number of rows of \emph{X}.}
  \item{p}{Number of columns of \emph{X}.}
}

\value{
Evaluate the trace (effective degree of freedom) of the product kernel smoother.}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\seealso{ \code{\link{ibr}}}


\keyword{smooth}
