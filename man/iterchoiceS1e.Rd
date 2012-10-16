\name{iterchoiceS1e}

\alias{iterchoiceS1e}

\title{Number of iterations selection for iterative bias reduction model}

\description{
Evaluate at each iteration proposed in the grid the value of different
criteria: GCV, AIC, corrected AIC, BIC and gMDL (along with the ddl and sigma
squared). The minimum of these
criteria gives an estimate of the optimal number of iterations.
This function is not intended to be used directly.
}

\usage{
iterchoiceS1e(y, K, tUy, eigenvaluesS1, ddlmini, ddlmaxi)}

\arguments{
  \item{y}{The response variable}
  \item{K}{A numeric vector which give the search grid for iterations}
  \item{eigenvaluesS1}{Vector of the eigenvalues of the
  symmetric smoothing matrix \emph{S}.}
  \item{tUy}{The transpose of the matrix of eigen vectors of the
  symmetric smoothing matrix \emph{S} times the vector of observation \emph{y}.}
  \item{ddlmini}{The number of eigen values of \emph{S} equal to 1.}
  \item{ddlmaxi}{The maximum df. No criteria are calculated beyond the
  number of iterations that leads to df bigger than this bound.}
}

\value{ Returns the values of GCV, AIC, corrected AIC, BIC, gMDL, df
and sigma squared for each value of the grid \code{K}. \code{Inf} are
returned if the iteration leads to a smoother with a df bigger than
\code{ddlmaxi}.}

\references{
  Cornillon, P.-A.; Hengartner, N.; Jegou, N. and Matzner-Lober, E. (2012)
  Iterative bias reduction: a comparative study.
\emph{Statistics and Computing}. Doi: 10.1007/s11222-012-9346-4
}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober}

\seealso{\code{\link{ibr}}, \code{\link{iterchoiceS1}}}


\keyword{smooth}
\keyword{multivariate}
