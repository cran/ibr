\name{iterchoiceA}

\alias{iterchoiceA}
\alias{critAgcv}
\alias{critAaic}
\alias{critAbic}
\alias{critAaicc}
\alias{critAgmdl}

\title{Selection of the number of iterations for iterative bias reduction smoothers}

\description{ The function \code{iterchoiceA} searches the interval from
\code{mini} to \code{maxi} for a minimum of the function
which calculates the chosen
\code{criterion} (\code{critAgcv}, \code{critAaic}, \code{critAbic},
\code{critAaicc} or \code{critAgmdl})  with respect to its first
argument (a given iteration \code{k}) using \code{\link{optimize}}. This function is not
intended to be used directly.  }

\usage{
iterchoiceA(n, mini, maxi, eigenvaluesA, tPADmdemiY, DdemiPA, 
ddlmini, ddlmaxi, y, criterion, fraction)}

\arguments{
  \item{n}{The number of observations.}
  \item{mini}{The lower end point of the interval to be searched.}
  \item{maxi}{The upper end point of the interval to be searched.}
  \item{eigenvaluesA}{Vector of the eigenvalues of the
  symmetric matrix \emph{A}.}
  \item{tPADmdemiY}{The transpose of the matrix of eigen vectors of the
  symmetric matrix \emph{A} times the inverse of the square root of the diagonal matrix \emph{D}.}
  \item{DdemiPA}{The square root of the diagonal matrix \emph{D} times
  the eigen vectors of the symmetric matrix \emph{A}.}
  \item{ddlmini}{The number of eigenvalues (numerically) equals to 1.}
  \item{ddlmaxi}{The  maximum df. No criterion is calculated and
  \code{Inf} is returned.}
  \item{y}{The vector of observations of dependant variable.}
  \item{criterion}{The criteria available are GCV (default, \code{"gcv"}),
  AIC (\code{"aic"}), corrected AIC  (\code{"aicc"}),   BIC
  (\code{"bic"}) or gMDL  (\code{"gmdl"}).}
  \item{fraction}{The subdivision of the interval [\code{mini},\code{maxi}].}
}
\details{See the reference for detailed explanation of \emph{A} and
  \emph{D}. The interval [\code{mini},\code{maxi}] is splitted into
  subintervals using \code{fraction}. In each subinterval the function
  \code{fcriterion} is minimzed  using \code{\link{optimize}} (with respect
     to its first argument) and the minimum (and its argument) of the
  result of these optimizations is returned.}

\value{
A list with components \code{iter} and \code{objective} which give the
(rounded) optimum number of iterations (between
\code{Kmin} and \code{Kmax}) and the value
     of the function at that real point (not rounded).}

\references{
  Cornillon, P.-A.; Hengartner, N.; Jegou, N. and Matzner-Lober, E. (2012)
  Iterative bias reduction: a comparative study.
\emph{Statistics and Computing}. Doi: 10.1007/s11222-012-9346-4
}
\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\seealso{\code{\link{ibr}}, \code{\link{iterchoiceA}}}


\keyword{smooth}
\keyword{multivariate}
