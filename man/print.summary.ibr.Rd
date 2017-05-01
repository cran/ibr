\name{print.summary.ibr}

\alias{print.summary.ibr}

\title{Printing iterative bias reduction summaries}

\description{
  \code{print} method for class  \dQuote{\code{summary.ibr}}.
}

\usage{
\method{print}{summary.ibr}(x,displaybw=FALSE, digits =
max(3, getOption("digits") - 3), \dots)}

\arguments{
  \item{x}{Object of class \code{\link{ibr}}.}
  \item{displaybw}{Boolean that indicates if bandwidth are printed or not.}
  \item{digits}{Rounds the values in its first argument to the specified
     number of significant digits.}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{
The function \code{print.summary.ibr} prints a list of summary
     statistics of the fitted  iterative bias reduction model given in \code{x}.}

\author{Pierre-Andre Cornillon, Nicolas Hengartner and Eric Matzner-Lober.}

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
\seealso{\code{\link{ibr}}, \code{\link{summary.ibr}}}

\examples{
\dontrun{data(ozone, package = "ibr")
res.ibr <- ibr(ozone[,-1],ozone[,1],df=1.2)
summary(res.ibr)
predict(res.ibr)}
}

\keyword{smooth}
\keyword{multivariate}
