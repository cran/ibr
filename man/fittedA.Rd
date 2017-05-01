\name{fittedA}

\alias{fittedA}

\title{Evaluates the fits for iterative bias reduction method}

\description{
Evaluates the fits for the iterative bias reduction smoother, using a kernel
smoother and its decomposition into a symmetric matrix and a diagonal
matrix. This function is not intended to be used directly.
}

\usage{
fittedA(n, eigenvaluesA, tPADmdemiY, DdemiPA, ddlmini, k)}

\arguments{
  \item{n}{The number of observations.}
  \item{eigenvaluesA}{Vector of the eigenvalues of the
  symmetric matrix \emph{A}.}
  \item{tPADmdemiY}{The transpose of the matrix of eigen vectors of the
  symmetric matrix \emph{A} times the inverse of the square root of the diagonal matrix \emph{D}.}
  \item{DdemiPA}{The square root of the diagonal matrix \emph{D} times
  the eigen vectors of the symmetric matrix \emph{A}.}
  \item{ddlmini}{The number of eigenvalues (numerically) equals to 1.}
  \item{k}{A scalar which gives the number of iterations.}
}
\details{See the reference for detailed explanation of \emph{A} and \emph{D}.}

\value{
Returns a list of two components: \code{fitted} contains fitted values
and \code{trace} contains the trace (effective degree of freedom) of the iterated
bias reduction smoother.}

\references{
  Cornillon, P.-A.; Hengartner, N.; Jegou, N. and Matzner-Lober, E. (2012)
  Iterative bias reduction: a comparative study.
\emph{Statistics and Computing}, \emph{23}, 777-791.

  Cornillon, P.-A.; Hengartner, N. and Matzner-Lober, E. (2013)
  Recursive bias estimation for multivariate regression smoothers Recursive 
bias estimation for multivariate regression smoothers.
\emph{ESAIM: Probability and Statistics}, \emph{18}, 483-502.

  Cornillon, P.-A.; Hengartner, N. and Matzner-Lober, E. (2017)
  Iterative Bias Reduction Multivariate Smoothing in R: The ibr Package.
\emph{Journal of Statistical Software}, \emph{77}, 1--26.
}

\author{Pierre-Andre Cornillon, Nicolas Hengartner and Eric Matzner-Lober.}

\seealso{\code{\link{ibr}}}


\keyword{smooth}
\keyword{multivariate}
