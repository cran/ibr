\name{ibr-package}
\alias{ibr-package}
\docType{package}
\title{
Iterative Bias Reduction
}
\description{
an R package for multivariate smoothing using Iterative Bias Reduction smoother.
}
\details{
  \itemize{\item We are interested in smoothing (the values of) a vector of
  \eqn{n} observations \eqn{y} by  \eqn{d} covariates measured at the
  same \eqn{n} observations (gathered in the matrix \eqn{X}).
 The iterated Bias Reduction produces a sequence of smoothers
  \deqn{\hat y=S_k y =(I - (I-S)^k)y,}
  where \eqn{S} is the pilot smoother which can be either a kernel or a
  thin plate spline smoother. In case of a kernel smoother, the kernel
  is built as a product of univariate kernels.
  \item The most important parameter of the iterated bias reduction is
  \eqn{k} the
  number of 
  iterationsr. Usually this parameter is unknown and is
  chosen from the search grid \code{K} to minimize the criterion (GCV,
  AIC, AICc, BIC or gMDL).\cr 
The user must choose the pilot smoother
  (kernel \code{"k"}, thin plate splines \code{"tps"} or Duchon splines \code{"ds"})
  plus the values of bandwidths (kernel)
  or \eqn{\lambda} thin plate splines). As the choice of these raw
  values depend on each particular dataset, one
  can rely on effective degrees of freedom or default values given as degree of
  freedom, see argument \code{df} of the main function \code{\link{ibr}}.}

\bold{Index of functions to be used by end user:}
\preformatted{
ibr:               Iterative bias reduction smoothing
plot.ibr:          Plot diagnostic for an ibr object
predict.ibr:       Predicted values using iterative bias reduction
                   smoothers
forward:           Variable selection for ibr (forward method)
print.summary.ibr: Printing iterative bias reduction summaries
summary.ibr:       Summarizing iterative bias reduction fits
}
}

  
\examples{
\dontrun{
data(ozone, package = "ibr")
res.ibr <- ibr(ozone[,-1],ozone[,1],smoother="k",df=1.1)
summary(res.ibr)
predict(res.ibr)
plot(res.ibr)
}
}


\author{
Pierre-Andre Cornillon, Nicolas Hengartner, Eric Matzner-Lober

Maintainer: Pierre-Andre Cornillon <pierre-andre.cornillon@supagro.inra.fr>
}
\keyword{ package }
