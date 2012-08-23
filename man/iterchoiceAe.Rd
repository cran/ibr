\name{iterchoiceAe}

\alias{iterchoiceAe}

\title{Selection of the number of iterations for iterative bias reduction smoothers}

\description{
Evaluates at each iteration proposed in the grid the value of different
criteria: GCV, AIC, corrected AIC, BIC and gMDL (along with the ddl and sigma
squared). The minimum of these
criteria gives an estimate of the optimal number of iterations.
This function is not intended to be used directly.
}

\usage{
iterchoiceAe(Y, K, eigenvaluesA, tPADmdemiY, DdemiPA, ddlmini,
ddlmaxi)}

\arguments{
  \item{Y}{The response variable.}
  \item{K}{A numeric vector which give the search grid for iterations.}
  \item{eigenvaluesA}{Vector of the eigenvalues of the
  symmetric matrix \emph{A}.}
  \item{tPADmdemiY}{The transpose of the matrix of eigen vectors of the
  symmetric matrix \emph{A} times the inverse of the square root of the diagonal matrix \emph{D}.}
  \item{DdemiPA}{The square root of the diagonal matrix \emph{D} times
  the eigen vectors of the symmetric matrix \emph{A}.}
  \item{ddlmini}{The number of eigenvalues (numerically) which are equal to 1.}
  \item{ddlmaxi}{The maximum df. No criteria are calculated beyond the
  number of iterations that leads to df bigger than this bound.}
}
\details{See the reference for detailed explanation of \emph{A} and \emph{D}}

\value{
Returns the values of GCV, AIC, corrected AIC, BIC, gMDL, df and sigma squared for each
value of the grid \code{K}. \code{Inf} are returned if the iteration leads
to a smoother with a df bigger than \code{ddlmaxi}.}

\references{
  Cornillon, P.-A.; Hengartner, N.; Jegou, N. and Matzner-Lober, E. (2012)
  Iterative bias reduction: a comparative study.
\emph{Statistics and Computing}. Doi: 10.1007/s11222-012-9346-4

  Cornillon, P.-A.; Hengartner, N. and Matzner-Lober, E. (2013)
  Recursive bias estimation for multivariate regression smoothers Recursive 
bias estimation for multivariate regression smoothers.
\emph{ESAIM: Probability and Statistics}. Doi: http://dx.doi.org/10.1051/ps/2013046
}
\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\seealso{\code{\link{ibr}}, \code{\link{iterchoiceA}}}


\keyword{smooth}
\keyword{multivariate}
