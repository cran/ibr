\name{print.summary.npregression}

\alias{print.summary.npregression}

\title{Printing iterative bias reduction summaries}

\description{
  \code{print} method for class  \dQuote{\code{summary.npregression}}.
}

\usage{
\method{print}{summary.npregression}(x,digits =
max(3, getOption("digits") - 3), \dots)}

\arguments{
  \item{x}{Object of class \code{\link{npregression}}.}
  \item{digits}{Rounds the values in its first argument to the specified
     number of significant digits.}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{
The function \code{print.summary.npregression} prints a list of summary
     statistics of the fitted  iterative bias reduction model given in \code{x}.}

\author{Pierre-Andre Cornillon, Nicolas Hengartner and Eric Matzner-Lober.}

\references{
Wand, M. P. and Jones, M. C. (1995). \emph{Kernel Smoothing}. Chapman and Hall, London.}

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
