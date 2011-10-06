\name{DuchonQ}

\alias{DuchonQ}

\title{Computes the semi-kernel of Duchon splines}

\description{
The function \code{DuchonQ} computes the semi-kernel of Duchon splines. This function is not intended to be used directly.
}

\usage{
DuchonQ(x,xk,m=2,s=0,symmetric=TRUE)}

\arguments{
  \item{x}{A numeric matrix of explanatory variables, with \emph{n} rows
  and \emph{p} columns.}
  \item{xk}{A numeric matrix of explanatory variables, with \emph{nk} rows
  and \emph{p} columns.}
  \item{m}{Order of derivatives.}
  \item{s}{Exponent for the weight function.}
  \item{symmetric}{Boolean: if \code{TRUE} only \code{x} is used and it
  computes the semi-kernel at observations of \code{x} (it should give the
  same result as \code{DuchonQ(x,xk,m,s,FALSE)}).}
}

\value{
The semi-kernel evaluated.}

\references{
Duchon, J. (1977) Splines minimizing rotation-invariant semi-norms in
Solobev spaces. in W. Shemp and K. Zeller (eds) \emph{Construction theory of
functions of several variables}, 85-100, Springer, Berlin. }

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\seealso{\code{\link{ibr}}}


\keyword{smooth}
\keyword{multivariate}
