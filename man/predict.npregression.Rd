\name{predict.npregression}

\alias{predict.npregression}

\title{Predicted values using using local polynomials

}

\description{
  Predicted values from a local polynomials of degree less than 2. See
   \code{\link[KernSmooth]{locpoly}} for fast binned implementation
  over an equally-spaced grid of local polynomial (gaussian kernel only) 
  \cr
 Missing values are not allowed.
}

\usage{
\method{predict}{npregression}(object, newdata, interval=
 c("none", "confidence", "prediction"), deriv=FALSE, \dots)}

\arguments{
  \item{object}{Object of class \code{\link{npregression}}.}
  \item{newdata}{An optional vector of values to be predicted. If omitted, the fitted values are used.}
  \item{interval}{Type of interval calculation. Only \code{none} is
  currently avalaible.}
  \item{deriv}{Bolean. If \code{TRUE} it returns the first derivative of
  the local  polynomial (of degree1).}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{ Produces a vector of predictions. If \code{deriv}
  is \code{TRUE} the value is a named list with components: \code{yhat}
  which contains predictions and (if relevant) \code{deriv} the
  first derivative of the local polynomial of degree 1.}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\references{
Wand, M. P. and Jones, M. C. (1995). \emph{Kernel Smoothing}. Chapman and Hall, London.}

\seealso{\code{\link{npregression}}, \code{\link{summary.npregression}},
  \code{\link[KernSmooth]{locpoly}}
}

\examples{
f <- function(x){sin(5*pi*x)}
n <- 100
x <- runif(n)
z <- f(x)
sigma2 <- 0.05*var(z)
erreur<-rnorm(n,0,sqrt(sigma2))
y<-z+erreur
grid <- seq(min(x),max(x),length=500)
res <- npregression(x,y,bandwidth=0.02,control.par=list(degree=1))
plot(x,y)
lines(grid,predict(res,grid))
}

\keyword{smooth}
\keyword{multivariate}
