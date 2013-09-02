\name{iterchoiceS1}

\alias{iterchoiceS1}
\alias{critS1gcv}
\alias{critS1aic}
\alias{critS1bic}
\alias{critS1aicc}
\alias{critS1gmdl}

\title{Number of iterations selection for iterative bias reduction model}

\description{
The function \code{iterchoiceS1} searches the interval from \code{mini} to
     \code{maxi} for a minimum of the function which calculates the chosen
\code{criterion} (\code{critS1gcv}, \code{critS1aic}, \code{critS1bic},
\code{critS1aicc} or \code{critS1gmdl})  with respect to its first
argument (a given iteration \code{k}) using \code{\link{optimize}}. This function is not intended to be used directly.
}

\usage{
iterchoiceS1(n, mini, maxi, tUy, eigenvaluesS1, ddlmini, ddlmaxi,
y, criterion, fraction)}

\arguments{
  \item{n}{The number of observations.}
  \item{mini}{The lower end point of the interval to be searched.}
  \item{maxi}{The upper end point of the interval to be searched.}
  \item{eigenvaluesS1}{Vector of the eigenvalues of the
  symmetric smoothing matrix \emph{S}.}
  \item{tUy}{The transpose of the matrix of eigen vectors of the
  symmetric smoothing matrix \emph{S} times the vector of observation \emph{y}.}
  \item{ddlmini}{The number of eigen values of \emph{S} equal to 1.}
  \item{ddlmaxi}{The  maximum df. No criterion is calculated and
  \code{Inf} is returned.}
  \item{y}{The vector of observations of dependant variable.}
  \item{criterion}{The criteria available are GCV (default, \code{"gcv"}),
  AIC (\code{"aic"}), corrected AIC  (\code{"aicc"}),   BIC
  (\code{"bic"}) or gMDL  (\code{"gmdl"}).}
  \item{fraction}{The subdivision of the interval [\code{mini},\code{maxi}].}
}

\details{The interval [\code{mini},\code{maxi}] is splitted into
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

  Cornillon, P.-A.; Hengartner, N. and Matzner-Lober, E. (2013)
  Recursive bias estimation for multivariate regression smoothers Recursive 
bias estimation for multivariate regression smoothers.
\emph{ESAIM: Probability and Statistics}. Doi: http://dx.doi.org/10.1051/ps/2013046
}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober}

\seealso{\code{\link{ibr}}, \code{\link{iterchoiceS1}}}


\keyword{smooth}
\keyword{multivariate}
