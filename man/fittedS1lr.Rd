\name{fittedS1lr}

\alias{fittedS1lr}

\title{Evaluate the fit for iterative bias reduction model}

\description{The function evaluates the fit for iterative bias reduction
  model for iteration \code{k}. This function is not intended to be used directly.}

\usage{
fittedS1lr(n,U,tUy,eigenvaluesS1,ddlmini,k,rank)
}

\arguments{
  \item{n}{The number of observations.}
  \item{U}{The the matrix of eigen vectors of the
  symmetric smoothing matrix \emph{S}.}
  \item{tUy}{The transpose of the matrix of eigen vectors of the
  symmetric smoothing matrix \emph{S} times the vector of observation \emph{y}.}
  \item{eigenvaluesS1}{Vector of the eigenvalues of the
  symmetric smoothing matrix \emph{S}.}
  \item{ddlmini}{The number of eigen values of \emph{S} equal to 1.}
  \item{k}{A numeric vector which gives the number of iterations}
 \item{rank}{The rank of lowrank splines.}
}
 
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

Wood, S.N. (2003) Thin plate regression
  splines. \emph{J. R. Statist. Soc. B}, \emph{65}, 95-114.
}
\details{see the reference for detailed explanation of computation of iterative bias reduction smoother}

\value{
Returns a vector containing the fit}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober}

\seealso{\code{\link{ibr}}}


\keyword{smooth}
