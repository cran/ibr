\name{poids}

\alias{poids}
\title{Product kernel evaluation}


\description{
Evaluate the product of kernel function at \emph{(X-valx)/bx}: Gaussian, Epanechnikov,
Uniform, Quartic. This function is not intended to be used directly.
}

\usage{
poids(kernelx,X,bx,valx,n,p)
}

\arguments{
  \item{kernelx}{Character string which allows to choose between gaussian kernel
  (\code{"g"}), Epanechnikov (\code{"e"}), uniform (\code{"u"}),
  quartic (\code{"q"}).}
  \item{X}{Matrix of explanatory variables, size \emph{n}, \emph{p}.}
  \item{bx}{The vector of different bandwidths, length \eqn{p}.}
  \item{valx}{The vector of length \eqn{p} at which the product kernel is evaluated.}
  \item{n}{Number of rows of \emph{X}.}
  \item{p}{Number of columns of \emph{X}.}
}

\value{
Returns a vector which coordinates are the values of
the product kernel at the given coordinate}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\seealso{ \code{\link{ibr}}}


\keyword{smooth}
