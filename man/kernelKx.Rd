\name{kernelKx}

\alias{kernelKx}

\title{Evaluates the smoothing matrix at x*}

\description{The function evaluates the kernel matrix 
  at arbitrary locations \emph{x*}.}

\usage{
kernelKx(X,Xetoile=NULL,bx,kernelx="g")
}

\arguments{
  \item{X}{Matrix of explanatory variables, size \emph{n}, \emph{p}.}
  \item{Xetoile}{Matrix of new design points \emph{x}* at which to predict the
  response variable, size \emph{n*}, \emph{p}. If \code{NULL} then
  \code{Xetoile} is equal to \code{X}.}
  \item{bx}{The vector of different bandwidths, length \eqn{p}.}
  \item{kernelx}{Character string which allows to choose between gaussian kernel
  (\code{"g"}), Epanechnikov (\code{"e"}), uniform (\code{"u"}),
  quartic (\code{"q"}).}
}

\value{
Returns the matrix \eqn{Kx} which \eqn{i,j} entry is equal to
\deqn{\prod_{k=1}^p K((x_{ik}^*-x_{jk})/h_k)} where \eqn{K(.)} is the
chosen kernel.
}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\seealso{\code{\link{kernelSx}}, \code{\link{ibr}}}


\keyword{smooth}
