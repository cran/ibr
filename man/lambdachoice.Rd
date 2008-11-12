\name{lambdachoice}

\alias{lambdachoice}

\title{Choice of bandwidth according to a given effective degree of freedom}

\description{
Perform a search for the different bandwidths in the given grid. For each explanatory
variable, the bandwidth is chosen such that the trace of the smoothing
matrix according to that variable (effective degree of freedom) is equal to a given value.
This function is not intended to be used directly.}

\usage{
lambdachoice(X,ddlobjectif,m=2,itermax)}

\arguments{
  \item{X}{A matrix with \eqn{n} rows (individuals) and \eqn{p} columns (numeric variables)}
  \item{ddlobjectif}{A numeric vector of length 1 which indicates the  desired effective degree of
  freedom (trace) of the smoothing   matrix for
    thin plate splines of order \code{m}.}
  \item{m}{Order of thin plate splines.}
  \item{itermax}{A scalar which controls the number of iterations for
    that search}
}

\value{
Returns the coefficient lambda that control smoothness for the desired
effective degree of freedom}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober}

\seealso{ \code{\link{ibr}}}


\keyword{smooth}
