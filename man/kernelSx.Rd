\name{kernelSx}

\alias{kernelSx}

\title{Evaluates the smoothing matrix at x*}

\description{The function evaluates the matrix of design weights to
  predict the response at arbitrary locations \emph{x}.
  This function is not intended to be used directly.}

\usage{
kernelSx(kernelx="g",X,Xetoile,bx)
}

\arguments{
  \item{kernelx}{Character string which allows to choose between gaussian kernel
  (\code{"g"}), Epanechnikov (\code{"e"}), uniform (\code{"u"}),
  quartic (\code{"q"}).}
  \item{X}{Matrix of explanatory variables, size \emph{n}, \emph{p}.}
  \item{Xetoile}{Matrix of new design points \emph{x}* at which to predict the
  response variable, size \emph{n*}, \emph{p}.}
  \item{bx}{The vector of different bandwidths, length \eqn{p}.}
}

\value{
Returns the matrix denoted in the paper by \eqn{Sx}, \emph{n*}, \emph{n}.}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\seealso{ \code{\link{ibr}}}


\keyword{smooth}
