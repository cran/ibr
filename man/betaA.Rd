\name{betaA}

\alias{betaA}

\title{Calculates coefficients for iterative bias reduction smoothers}

\description{
Calculates the coefficients for the iterative bias reduction  smoothers. This function is not intended to be used directly.}

\usage{
betaA(n, eigenvaluesA, tPADmdemiY, DdemiPA, ddlmini, k, index0)}

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
  \item{index0}{The index of the first eigen values of \emph{S}
  numerically equal to 0.}
}
\details{See the reference for detailed explanation of \emph{A} and
  \emph{D} and the meaning of coefficients.}

\value{
Returns the  vector of coefficients (of length \emph{n}, the number of observations.)}

\references{
  Cornillon, P.-A.; Hengartner, N.; Jegou, N. and Matzner-Lober, E. (2012)
  Iterative bias reduction: a comparative study.
\emph{Statistics and Computing}, \emph{23}, 777-791.

  Cornillon, P.-A.; Hengartner, N. and Matzner-Lober, E. (2013)
  Recursive bias estimation for multivariate regression smoothers Recursive 
bias estimation for multivariate regression smoothers.
\emph{ESAIM: Probability and Statistics}, \emph{18}, 483-502.
}
\author{Pierre-Andre Cornillon, Nicolas Hengartner and Eric Matzner-Lober.}

\seealso{\code{\link{ibr}}}


\keyword{smooth}
\keyword{multivariate}
