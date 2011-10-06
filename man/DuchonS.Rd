\name{DuchonS}

\alias{DuchonS}

\title{Computes the semi-kernel of Duchon splines}

\description{
The function \code{DuchonS} computes the semi-kernel of Duchon splines. This function is not intended to be used directly.
}

\usage{
DuchonS(x,m=2)}

\arguments{
  \item{x}{A numeric matrix of explanatory variables, with \emph{n} rows
  and \emph{p} columns.}
  \item{m}{Order of derivatives.}
}

\value{
The polynomial part evaluated.}

\references{
Duchon, J. (1977) Splines minimizing rotation-invariant semi-norms in
Solobev spaces. in W. Shemp and K. Zeller (eds) \emph{Construction theory of
functions of several variables}, 85-100, Springer, Berlin. }

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\seealso{\code{\link{ibr}}}


\keyword{smooth}
\keyword{multivariate}
