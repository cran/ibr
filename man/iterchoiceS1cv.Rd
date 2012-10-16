\name{iterchoiceS1cv}

\alias{iterchoiceS1cv}

\title{Selection of the number of iterations for iterative bias
  reduction smoothers with base thin-plate splines or duchon splines smoother}

\description{
The function \code{iterchoiceS1cv} searches the interval from \code{mini} to
     \code{maxi} for a minimum of the function \code{criterion} with respect
     to its first argument using \code{\link{optimize}}. This function is not intended to be used directly.
}

\usage{
iterchoiceS1cv(X, y, lambda, df, ddlmini, ntest, ntrain,
Kfold, type, npermut, seed, Kmin, Kmax, criterion, m, s,
fraction)}

\arguments{
  \item{X}{A numeric matrix of explanatory variables, with \emph{n} rows
  and \emph{p} columns.}
  \item{y}{A numeric vector of variable to be explained of length \emph{n}.}
  \item{lambda}{A numeric positive coefficient that governs the
  amount of penalty (coefficient lambda).}
  \item{df}{A numeric vector of length 1 which is multiplied by the minimum df of thin
  plate splines ; This argument is useless if
  \code{lambda} is supplied (non null).}
  \item{ddlmini}{The number of eigenvalues equals to 1.}
  \item{ntest}{The number of observations in test set.}
  \item{ntrain}{The number of observations in training set.}
  \item{Kfold}{Either the number of folds or a boolean or \code{NULL}.}
  \item{type}{A character string in
  \code{random},\code{timeseries},\code{consecutive}, \code{interleaved}
  and give the type of segments.}
  \item{npermut}{The number of random draw (with replacement), used for
    \code{type="random"}.}
  \item{seed}{Controls the seed of random generator
    (via \code{\link{set.seed}}).}
  \item{Kmin}{The minimum number of bias correction  iterations of the
  search grid considered by
  the model selection procedure for selecting the optimal number of iterations.}
  \item{Kmax}{The maximum number of bias correction  iterations of the
  search grid considered by
  the model selection procedure for selecting the optimal number of iterations.}
  \item{criterion}{The criteria available are map (\code{"map"}) or rmse
    (\code{"rmse"}).}
  \item{m}{The order of derivatives for the penalty (for thin plate
  splines it is the order). This integer \emph{m} must verify
  2\emph{m}+2\emph{s}/\emph{d}>1, where \emph{d} is the number of explanatory
  variables.} 
  \item{s}{The power of weighting function. For thin plate splines
  \emph{s} is equal to 0. This real must be strictly smaller than \emph{d}/2
  (where \emph{d} is the number of explanatory  variables) and must
  verify 2\emph{m}+2\emph{s}/\emph{d}. To get pseudo-cubic splines,
  choose \emph{m}=2 and \emph{s}=(\emph{d}-1)/2 (See Duchon, 1977).} 
  \item{fraction}{The subdivision of the interval [\code{Kmin},\code{Kmax}].}
}

\value{
Returns the optimum number of iterations (between \code{Kmin} and \code{Kmax}).}

\references{
  Cornillon, P.-A.; Hengartner, N.; Jegou, N. and Matzner-Lober, E. (2012)
  Iterative bias reduction: a comparative study.
\emph{Statistics and Computing}. Doi: 10.1007/s11222-012-9346-4

Duchon, J. (1977) Splines minimizing rotation-invariant semi-norms in
Solobev spaces. in W. Shemp and K. Zeller (eds) \emph{Construction theory of
functions of several variables}, 85-100, Springer, Berlin. }

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\seealso{\code{\link{ibr}}}


\keyword{smooth}
\keyword{multivariate}
