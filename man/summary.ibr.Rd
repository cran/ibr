\name{summary.ibr}

\alias{summary.ibr}

\title{Summarizing iterative bias reduction fits}

\description{
  \code{summary} method for class  \dQuote{\code{ibr}}.
}

\usage{
\method{summary}{ibr}(object,  criteria="call", \dots)}

\arguments{
  \item{object}{Object of class \code{\link{ibr}}.}
  \item{criteria}{Character string which gives the criteria evaluated for the model. The criteria available are GCV (default, \code{"gcv"}),
  AIC (\code{"aic"}), corrected AIC  (\code{"aicc"}),   BIC
  (\code{"bic"}) or gMDL  (\code{"gmdl"}). The string \code{"call"}
  return the criterion used in the call of \code{ibr}.}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{
The function \code{summary.ibr} computes and returns a list of summary
     statistics of the fitted  iterative bias reduction smoother given in \code{object}}

\author{Pierre-Andre Cornillon, Nicolas Hengartner and Eric Matzner-Lober.}

\references{
Cornillon, P. A., Hengartner, N. and Matzner-Lober,  E. (2009) Recursive
Bias Estimation for high dimensional regression smoothers. \emph{submitted}.
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