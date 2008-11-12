\name{calcA}

\alias{calcA}

\title{Decomposition of the kernel smoother}

\description{
Calculates the decomposition of the kernel smoothing matrix in two part: a
diagonal matrix \emph{D} and a symmetric matrix \emph{A}.
This function is not intended to be used directly.
}

\usage{
calcA(X,bx,kernelx="g")}

\arguments{
  \item{X}{The matrix of explanatory variables, size \emph{n}, \emph{p}.}
  \item{bx}{The vector of bandwidth of length \emph{p}.}
  \item{kernelx}{Character string which allows to choose between gaussian kernel
  (\code{"g"}), Epanechnikov (\code{"e"}), uniform (\code{"u"}),
  quartic (\code{"q"}).}
}
\details{see the reference for detailed explanation of \emph{A} and
  \emph{D} and the meaning of coefficients.}

\value{
  Returns a list containing two matrices: the symmetric matrix \emph{A}
  in component \code{A}) and the square root of the diagonal matrix
  \emph{D} in the component \code{Ddemi} and the trace of the smoother
  in the component \code{df}.}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\seealso{\code{\link{ibr}}}


\keyword{smooth}
\keyword{multivariate}
