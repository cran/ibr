\name{iterchoiceS1lrcve}

\alias{iterchoiceS1lrcve}

\title{Selection of the number of iterations for iterative bias
  reduction smoothers with base lowrank thin-plate splines smoother or duchon splines smoother}

\description{
Evaluates at each iteration proposed in the grid the cross-validated
root mean squared error (RMSE) and mean of the relative absolute error (MAP). The minimum of these
criteria gives an estimate of the optimal number of iterations.
This function is not intended to be used directly.
}

\usage{
iterchoiceS1lrcve(X, y, lambda, rank, bs, listvarx, df, ddlmini, ntest, ntrain,
Kfold, type, npermut, seed, Kmin, Kmax, m, s)}

\arguments{
  \item{X}{A numeric matrix of explanatory variables, with \emph{n} rows
  and \emph{p} columns.}
  \item{y}{A numeric vector of variable to be explained of length \emph{n}.}
  \item{lambda}{A numeric positive coefficient that governs the
  amount of penalty (coefficient lambda).}
 \item{rank}{The rank of lowrank splines.}
 \item{bs}{The type rank of lowrank splines: \code{tps} or \code{ds}.}
 \item{listvarx}{The vector of the names of explanatory variables}
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
  \item{m}{The order of derivatives for the penalty (for thin plate
  splines it is the order). This integer \emph{m} must verify
  2\emph{m}+2\emph{s}/\emph{d}>1, where \emph{d} is the number of explanatory
  variables.} 
  \item{s}{The power of weighting function. For thin plate splines
  \emph{s} is equal to 0. This real must be strictly smaller than \emph{d}/2
  (where \emph{d} is the number of explanatory  variables) and must
  verify 2\emph{m}+2\emph{s}/\emph{d}. To get pseudo-cubic splines,
  choose \emph{m}=2 and \emph{s}=(\emph{d}-1)/2 (See Duchon).} 
}

\value{
Returns the values of RMSE and MAP for each
value of the grid \code{K}. \code{Inf} are returned if the iteration leads
to a smoother with a df bigger than \code{ddlmaxi}.}

\references{   Cornillon, P.-A.; Hengartner, N.; Jegou, N. and Matzner-Lober, E. (2012)
  Iterative bias reduction: a comparative study.
\emph{Statistics and Computing}, \emph{23}, 777-791.



  Cornillon, P.-A.; Hengartner, N. and Matzner-Lober, E. (2013)
  Recursive bias estimation for multivariate regression smoothers Recursive 
bias estimation for multivariate regression smoothers.
\emph{ESAIM: Probability and Statistics}, \emph{18}, 483-502.

  Cornillon, P.-A.; Hengartner, N. and Matzner-Lober, E. (2017)
  Iterative Bias Reduction Multivariate Smoothing in R: The ibr Package.
\emph{Journal of Statistical Software}, \emph{77}, 1--26.

Duchon, J. (1977) Splines minimizing rotation-invariant semi-norms in
Solobev spaces. in W. Shemp and K. Zeller (eds) \emph{Construction theory of
  functions of several variables}, 85-100, Springer, Berlin.

Wood, S.N. (2003) Thin plate regression
  splines. \emph{J. R. Statist. Soc. B}, \emph{65}, 95-114.
}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\seealso{\code{\link{ibr}}}


\keyword{smooth}
\keyword{multivariate}
