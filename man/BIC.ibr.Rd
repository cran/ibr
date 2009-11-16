\name{BIC}
\alias{BIC}
\alias{GCV}
\alias{AICc}
\alias{BIC.ibr}
\alias{GCV.ibr}
\alias{AICc.ibr}

\title{Information Criterion for ibr}

\description{
  Functions calculating the  Bayesian Informative Criterion , the
  Generalized Cross Validation criterion and the Corrected Akaike information criterion.
}

\usage{
\method{BIC}{ibr}(object, \dots)

\method{GCV}{ibr}(object, \dots)

\method{AICc}{ibr}(object, \dots)}

\arguments{
  \item{object}{A fitted model object of class ibr.}
  \item{\dots}{Only for compatibility purpose with \code{BIC} of
  \code{nlme} package.}
}
\details{
  The ibr method for \code{BIC}, \code{BIC.ibr()} calculates
  \eqn{\log(sigma^2)+log(n)*df/n}{log(sigma^2)+log(n)*df/n}, where \emph{df} is the trace
  of the smoother.

The ibr method for \code{GCV}, \code{GCV.ibr()} calculates
  \eqn{\log(sigma^2)-2*\log(1-df/n)}{log(sigma^2)-log(1-*df/n)}
  
The ibr method for \code{AICc}, \code{AICc.ibr()} calculates
  \eqn{\log(sigma^2)+1+(2*(df+1))/(n-df-2)}{log(sigma^2)+1+(2*(df+1))/(n-df-2)}.
  
}
\value{
  Returns a numeric value
  with the corresponding BIC, GCV or AICc.
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
res.ibr <- ibr(ozone[,-1],ozone[,1])
BIC(res.ibr)
GCV(res.ibr)
AICc(res.ibr)
}
}

\keyword{smooth}
\keyword{multivariate}
