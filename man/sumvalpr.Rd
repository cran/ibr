\name{sumvalpr}

\alias{sumvalpr}

\title{Sum of a geometric series}

\description{
Calculates the sum of the first (k+1) terms of a geometric series with initial term 1
  and common ratio equal to \code{valpr} (lower or equal to 1).}

\usage{
sumvalpr(k,n,valpr,index1,index0)}

\arguments{
  \item{k}{The number of terms minus 1.}
  \item{n}{The length of \code{valpr}.}
  \item{valpr}{Vector of common ratio in decreasing order.}
  \item{index1}{The index of the last common ratio equal to 1.}
  \item{index0}{The index of the first common ratio equal to 0.}
}

\value{
Returns the vector of the sums of the first (k+1) terms of the geometric series.}

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
