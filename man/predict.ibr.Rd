\name{predict.ibr}

\alias{predict.ibr}

\title{Predicted values using iterative bias reduction smoothers}

\description{
  Predicted values from iterative bias reduction object.
  \cr
 Missing values are not allowed.
}

\usage{
\method{predict}{ibr}(object, newdata, interval=
 c("none", "confidence", "prediction"), \dots)}

\arguments{
  \item{object}{Object of class \code{\link{ibr}}.}
  \item{newdata}{An optional matrix in which to look for variables with which to predict. If omitted, the fitted values are used.}
  \item{interval}{Type of interval calculation. Only \code{none} is
  currently avalaible.}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{
  Produces a vector of predictions.}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\references{
Cornillon, P. A., Hengartner, N. and Matzner-Lober,  E. (2009) Recursive
Bias Estimation for high dimensional regression
smoothers. \emph{submitted}.
}

\seealso{\code{\link{ibr}}, \code{\link{summary.ibr}}}

\examples{
\dontrun{data(ozone, package = "ibr")
res.ibr <- ibr(ozone[,-1],ozone[,1],df=1.2,K=1:500)
summary(res.ibr)
predict(res.ibr)}
}

\keyword{smooth}
\keyword{multivariate}
