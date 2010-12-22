\name{summary.npregression}

\alias{summary.npregression}

\title{Summarizing local polynomial fits}

\description{
  \code{summary} method for class  \dQuote{\code{npregression}}.
}

\usage{
\method{summary}{npregression}(object,  criteria="call", \dots)}

\arguments{
  \item{object}{Object of class \code{\link{npregression}}.}
  \item{criteria}{Character string which gives the criteria evaluated for the model. The criteria available are GCV (default, \code{"gcv"}),
  AIC (\code{"aic"}), corrected AIC  (\code{"aicc"}),   BIC
  (\code{"bic"}) or gMDL  (\code{"gmdl"}). The string \code{"call"}
  return the criterion used in the call of \code{npregression}.}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{
The function \code{summary.npregression} computes and returns a list of summary
     statistics of the local polynomial smoother given in \code{object}}

\author{Pierre-Andre Cornillon, Nicolas Hengartner and Eric Matzner-Lober.}

\references{
Wand, M. P. and Jones, M. C. (1995). \emph{Kernel Smoothing}. Chapman and Hall, London.}

\seealso{\code{\link{npregression}}, \code{\link{summary.npregression}}}

\examples{
f <- function(x){sin(5*pi*x)}
n <- 100
x <- runif(n)
z <- f(x)
sigma2 <- 0.05*var(z)
erreur <- rnorm(n,0,sqrt(sigma2))
y <- z+erreur
res <- npregression(x,y,bandwidth=0.02)
summary(res)
}

\keyword{smooth}
\keyword{multivariate}
