\name{lambdachoicelr}

\alias{lambdachoicelr}

\title{Choice of bandwidth according to a given effective degree of freedom}

\description{
Perform a search for the different bandwidths in the given grid. For each explanatory
variable, the bandwidth is chosen such that the trace of the smoothing
matrix according to that variable (effective degree of freedom) is equal to a given value.
This function is not intended to be used directly.}

\usage{
lambdachoicelr(x,ddlobjectif,m=2,s=0,rank,itermax,bs,listvarx)}

\arguments{
  \item{x}{A matrix with \eqn{n} rows (individuals) and \eqn{p} columns (numeric variables)}
  \item{ddlobjectif}{A numeric vector of length 1 which indicates the  desired effective degree of
  freedom (trace) of the smoothing   matrix for
    thin plate splines of order \code{m}.}
  \item{m}{The order of derivatives for the penalty (for thin plate
  splines it is the order). This integer \emph{m} must verify
  2\emph{m}+2\emph{s}/\emph{d}>1, where \emph{d} is the number of explanatory
  variables.} 
  \item{s}{The power of weighting function. For thin plate splines
  \emph{s} is equal to 0. This real must be strictly smaller than \emph{d}/2
  (where \emph{d} is the number of explanatory  variables) and must
  verify 2\emph{m}+2\emph{s}/\emph{d}. To get pseudo-cubic splines,
  choose \emph{m}=2 and \emph{s}=(\emph{d}-1)/2 (See Duchon, 1977).}
  \item{itermax}{A scalar which controls the number of iterations for
    that search}
 \item{rank}{The rank of lowrank splines.}
 \item{bs}{The type rank of lowrank splines: \code{tps} or \code{ds}.}
 \item{listvarx}{The vector of the names of explanatory variables}
}

\references{
Duchon, J. (1977) Splines minimizing rotation-invariant semi-norms in
Solobev spaces. in W. Shemp and K. Zeller (eds) \emph{Construction theory of
functions of several variables}, 85-100, Springer, Berlin. 

Wood, S.N. (2003) Thin plate regression
  splines. \emph{J. R. Statist. Soc. B}, \emph{65}, 95-114.
}
\value{
Returns the coefficient lambda that control smoothness for the desired
effective degree of freedom}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober}

\seealso{ \code{\link{ibr}}}


\keyword{smooth}
