\name{bwchoice}

\alias{bwchoice}

\title{Choice of bandwidth achieving a prescribed effective degree of freedom}

\description{
Perform a search for the bandwidths in the given grid. For each explanatory
variable, the bandwidth is chosen such that the trace of the smoothing
matrix according to that variable (effective degree of freedom) is equal to a
prescribed value.
This function is not intended to be used directly.}

\usage{
bwchoice(X,objectif,kernelx="g",itermax=1000)}

\arguments{
  \item{X}{A matrix with \eqn{n} rows (individuals) and \eqn{p} columns (numeric variables).}
  \item{objectif}{A numeric vector of either length 1 or length equal to the
  number of columns of \code{X}. It indicates the  desired effective degree of
  freedom (trace) of the smoothing   matrix for
    each variable. \code{objectif} is repeated when the length of vector
  \code{objectif} is 1.}
  \item{kernelx}{String which allows to choose between gaussian kernel
  (\code{"g"}), Epanechnikov (\code{"e"}), uniform (\code{"u"}),
  quartic (\code{"q"}).}
  \item{itermax}{A scalar which controls the number of iterations for
    that search.}
}

\value{
Returns a vector of length d, the number of explanatory variable, where
each coordinate is the value of the selected bandwidth for each
explanatory variable}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\seealso{ \code{\link{ibr}}}


\keyword{smooth}
