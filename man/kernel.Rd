\name{kernel}

\alias{kernel}
\alias{gaussien}
\alias{epane}
\alias{uniform}
\alias{quartic}
\title{Kernel evaluation}


\description{
Evaluate the kernel function at \emph{x}: Gaussian, Epanechnikov,
Uniform, Quartic. This function is not intended to be used directly.
}

\usage{
gaussien(X)
epane(X)
uniform(X)
quartic(X)
}

\arguments{
  \item{X}{The value where the function has to be evaluate, should be a
    numeric and can be a scalar, a vector or a matrix}
}

\value{
Returns a scalar, a vector or a matrix which coordinates are the values of
the kernel at the given coordinate}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\seealso{ \code{\link{ibr}}}


\keyword{smooth}
