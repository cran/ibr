\name{plot.forwardibr}

\alias{plot.forwardibr}

\title{Plot diagnostic for an ibr object}

\description{
  One plot  is  currently available: a plot
     of residuals against fitted values.
}

\usage{
\method{plot}{forwardibr}(x,\ldots )} 

\arguments{
  \item{x}{Object of class \code{\link{forwardibr}}.}
  \item{\dots}{further arguments passed to \code{\link{image}}.}
}

\value{ The function \code{plot.forwardibr} give an image plot of the
values of the criterion obtained by the forward selection process. Image
is read from the bottom to the top. At the bottom row, there are all the
univariate models and the selected variable is given by the lowest
criterion. This variable is selected for the second row. At the second
(bottom) row the second variable included is those which give the lowest
criterion for this row etc. All the variables included in the final
model (selected by forward search) are numbered on the image (by order of
inclusion).}

\author{Pierre-Andre Cornillon, Nicolas Hengartner and Eric Matzner-Lober.}

\references{
Cornillon, P. A., Hengartner, N. and Matzner-Lober,  E. (2009) Recursive
Bias Estimation for high dimensional regression smoothers. \emph{submitted}.
}
\seealso{\code{\link{ibr}}, \code{\link{forward}}}

\examples{
\dontrun{data(ozone, package = "ibr")
ibrsel <- forward(ibr(ozone[,-1],ozone[,1],df=1.2)
plot(ibrsel)
plot(apply(ibrsel,1,min,na.rm=TRUE),type="l")
}
}

\keyword{smooth}
\keyword{multivariate}
