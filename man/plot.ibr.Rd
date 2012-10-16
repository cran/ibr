\name{plot.ibr}

\alias{plot.ibr}

\title{Plot diagnostic for an ibr object}

\description{
  One plot  is  currently available: a plot
     of residuals against fitted values.
}

\usage{
\method{plot}{ibr}(x,\ldots )} 

\arguments{
  \item{x}{Object of class \code{\link{ibr}}.}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{
The function \code{plot.ibr} computes and returns a list of summary
     statistics of the fitted  iterative bias reduction smoother given in \code{object}}

\author{Pierre-Andre Cornillon, Nicolas Hengartner and Eric Matzner-Lober.}

\references{
  Cornillon, P.-A.; Hengartner, N.; Jegou, N. and Matzner-Lober, E. (2012)
  Iterative bias reduction: a comparative study.
\emph{Statistics and Computing}. Doi: 10.1007/s11222-012-9346-4
}
\seealso{\code{\link{ibr}}, \code{\link{summary.ibr}}}

\examples{
\dontrun{data(ozone, package = "ibr")
res.ibr <- ibr(ozone[,-1],ozone[,1],df=1.2)
plot(res.ibr)}
}

\keyword{smooth}
\keyword{multivariate}
