\name{npregress}

\alias{npregress}
\alias{print.npregress}
\alias{residuals.npregress}

\title{Local polynomials smoothing}

\description{
  Predicted values from a local polynomials of degree less than 2.  \cr
Missing values are not allowed.
}

\usage{
npregress(x, y, criterion="rmse", bandwidth=NULL,kernel="g",
             control.par=list(), cv.options=list())}

\arguments{
  \item{x}{A numeric vector of explanatory variable of length \emph{n}.}
  \item{y}{A numeric vector of variable to be explained of length \emph{n}.}
  \item{criterion}{Character string. If the bandwidth 
    (\code{bandwidth}) is missing or 
    \code{NULL} the number of iterations is chosen using
    \code{criterion}. The criterion available is (cross-validated) rmse
  (\code{"rmse"}) and mean (relative) absolute error.}
\item{bandwidth}{The kernel bandwidth smoothing parameter (a numeric vector of either length 1).}
\item{kernel}{Character string which allows to choose between gaussian kernel
  (\code{"g"}), Epanechnikov (\code{"e"}), uniform (\code{"u"}),
  quartic (\code{"q"}).}
\item{control.par}{A named list that control optional parameters. The
  two components are \code{bandwidth} for compatibility with \code{\link{ibr}}
  arguments and \code{degree} which controls the degree of the local
  polynomial regression. If argument \code{bandwidth} is not null or missing, its
  value is used instead \code{control.par$bandwidth}. \code{degree} must
  be smaller than 2. For (gaussian binned) local polynomial see
  \code{\link[KernSmooth]{locpoly}}}
  \item{cv.options}{A named list which controls the way to do cross
  validation with component \code{gridbw},
  \code{ntest}, \code{ntrain}, \code{Kfold}, \code{type},
  \code{seed}, \code{method} and \code{npermut}. \code{gridbw} is
  numeric vector which contains the search grid for optimal bandwidth  (default
  to \code{1/n*(1+1/n)^(0:kmax)}, with \code{kmax=floor(log(n*diff(range(x))/3)/log(1+1/n))}). \code{ntest} is the number of observations in test set and
  \code{ntrain} is the number of observations in training set. Actually,
  only one of these is needed the other can be \code{NULL} or missing. \code{Kfold} a boolean or an integer. If
  \code{Kfold} is \code{TRUE} then the number of fold is deduced from
  \code{ntest} (or \code{ntrain}).  \code{type} is a character string in
  \code{random},\code{timeseries},\code{consecutive}, \code{interleaved}
  and give the type of segments.  \code{seed} controls the seed of
  random generator. \code{npermut} is the number of random draws.   If
  \code{cv.options} is \code{list()}, then component \code{ntest} is set to
  \code{1}, \code{type} is consecutive, \code{Kfold} is \code{TRUE}, and
  the other components are \code{NULL}, which leads to leave-one-out
  cross-validation.}
}

\value{
Returns an object of class \code{npregress} which is a list including:
  \item{bandwidth}{The kernel bandwidth smoothing parameter.}
  \item{residuals}{Vector of residuals.}
  \item{fitted}{Vector of fitted values.}
  \item{df}{The effective degree of freedom of the smoother.}
   \item{call}{A list containing four components: \code{x} contains the
    initial explanatory variables, \code{y} contains the
    initial dependant variables, 
    \code{criterion} contains the chosen criterion, \code{kernel} the
    kernel and \code{degree} the chosen degree}
  \item{criteria}{either a named list containing the bandwidth search
    grid and all the criteria (\code{rmse} and \code{mae}) evaluated on the
    grid \code{gridbw}. If the bandwidth
    \code{bandwidth} is given by the user \code{NULL} is returned}
}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\note{See \code{\link[KernSmooth]{locpoly}} for fast binned implementation
  over an equally-spaced grid of local polynomial. See \code{\link{ibr}}
  for univariate and \strong{multivariate} smoothing. 
}
\references{
Wand, M. P. and Jones, M. C. (1995). \emph{Kernel Smoothing}. Chapman and Hall, London.}

\seealso{\code{\link{predict.npregress}},
  \code{\link{summary.npregress}},
  \code{\link[KernSmooth]{locpoly}}, \code{\link{ibr}}}

\examples{
f <- function(x){sin(5*pi*x)}
n <- 100
x <- runif(n)
z <- f(x)
sigma2 <- 0.05*var(z)
erreur <- rnorm(n,0,sqrt(sigma2))
y <- z+erreur
res <- npregress(x,y,bandwidth=0.02)
summary(res)
ord <- order(x)
plot(x,y)
lines(x[ord],predict(res)[ord])
}

\keyword{smooth}
\keyword{multivariate}
