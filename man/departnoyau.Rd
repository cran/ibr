\name{departnoyau}

\alias{departnoyau}
\title{Trace of the product kernel smoother}


\description{
 Search bandwidth for each univariate kernel smoother such that the
 product of these univariate kernel gives a kernel smoother with a
 chosen effective degree of freedom (trace of the smoother). The bandwidths are
 constrained to give, for each explanatory
variable, a kernel smoother with same trace as the others.
This function is not intended to be used directly.
}

\usage{
departnoyau(df, x, kernel, dftobwitmax, n, p, dfobjectif)
}

\arguments{
  \item{df}{A numeric vector giving the effective degree of freedom (trace) of the
    univariate smoothing  matrix  for each variable of 
    \code{x}.}
  \item{x}{Matrix of explanatory variables, size \emph{n}, \emph{p}.}
  \item{kernel}{Character string which allows to choose between gaussian kernel
    (\code{"g"}), Epanechnikov (\code{"e"}), uniform (\code{"u"}),
    quartic (\code{"q"}).}
  \item{dftobwitmax}{Specifies the maximum number of iterations
    transmitted   to \code{\link{uniroot}} function.}
  \item{n}{Number of rows of data matrix \emph{x}.}
  \item{p}{Number of columns of data matrix \emph{x}.}
  \item{dfobjectif}{A numeric vector of length 1 which indicates
    the  desired effective degree of freedom (trace) of the smoothing
    matrix (product kernel smoother) for \code{x}.}
 }

\value{
Returns the desired bandwidths.}


\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\seealso{ \code{\link{ibr}}}


\keyword{smooth}
