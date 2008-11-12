\name{AIC.ibr}

\alias{AIC.ibr}

\title{Summarizing iterative bias reduction fits}

\description{
  Generic function calculating the Akaike information criterion for
  one model objects of ibr class for which a log-likelihood value
  can be obtained, according to the formula
  \eqn{-2 \log(sigma^2) + k df/n}{-2 log(sigma^2) + k df/n},
  where \eqn{df}{df} represents the effective degree of freedom (trace) of the
  smoother in the
  fitted model, and \eqn{k = 2} for the usual AIC, or \eqn{k = \log(n)}
  (\eqn{n} the number of observations) for the so-called BIC or SBC
  (Schwarz's Bayesian criterion).
}

\usage{
\method{AIC}{ibr}(object, \dots, k = 2)}

\arguments{
  \item{object}{A fitted model object of class ibr.}
  \item{\dots}{Not used.}
  \item{k}{Numeric, the \emph{penalty} per parameter to be used; the
    default \code{k = 2} is the classical AIC.}
}
\details{
  The ibr method for \code{AIC}, \code{AIC.ibr()} calculates
  \eqn{\log(sigma^2)+2*df/n}{log(sigma^2)+2*df/n}, where \emph{df} is the trace
  of the smoother.
}
\value{
  returns a numeric value
  with the corresponding AIC (or BIC, or \dots, depending on \code{k}).
}
\author{Pierre-Andre Cornillon, Nicolas Hengartner and Eric Matzner-Lober.}

\references{
Hurvich, C. M., Simonoff J. S. and Tsai, C. L. (1998) Smoothing
Parameter Selection in Nonparametric Regression Using an Improved Akaike
Information Criterion. \emph{Journal of the Royal Statistical Society, Series B}, 60, 271-293 . 
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
