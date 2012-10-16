\name{plot.forwardibr}

\alias{plot.forwardibr}

\title{Plot diagnostic for an ibr object}

\description{
  One plot  is  currently available: a plot
     of residuals against fitted values.
}

\usage{
\method{plot}{forwardibr}(x,global=FALSE,\ldots )} 

\arguments{
  \item{x}{Object of class \code{\link{forwardibr}}.}
  \item{global}{Boolean: if \code{global} is \code{TRUE} the color code
  is between the min and the max of \code{x} (except infinite value); if
  \code{global} is \code{FALSE} the color code is between the min and
  the max of each row.}
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
  Cornillon, P.-A.; Hengartner, N.; Jegou, N. and Matzner-Lober, E. (2012)
  Iterative bias reduction: a comparative study.
\emph{Statistics and Computing}. Doi: 10.1007/s11222-012-9346-4
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
