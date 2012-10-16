\name{iterchoiceAcve}

\alias{iterchoiceAcve}

\title{Selection of the number of iterations for iterative bias reduction smoothers}

\description{
Evaluates at each iteration proposed in the grid the cross-validated
root mean squared error (RMSE) and mean of the relative absolute error (MAP). The minimum of these
criteria gives an estimate of the optimal number of iterations.
This function is not intended to be used directly.
}

\usage{
iterchoiceAcve(X, y, bx, df, kernelx, ddlmini, ntest, ntrain,
Kfold, type, npermut, seed, Kmin, Kmax)}

\arguments{
  \item{X}{A numeric matrix of explanatory variables, with \emph{n} rows
  and \emph{p} columns.}
  \item{y}{A numeric vector of variable to be explained of length \emph{n}.}
  \item{bx}{The vector of different bandwidths, length \eqn{p}.}
  \item{df}{A numeric vector of either length 1 or length equal to the
    number of columns of \code{x}. If \code{smoother="k"}, it indicates
    the  desired effective degree of
  freedom (trace) of the smoothing   matrix for
    each variable ; \code{df} is repeated when the length of vector
  \code{df} is 1. This argument is useless if
  \code{bandwidth} is supplied (non null).}
  \item{kernelx}{Character string which allows to choose between gaussian kernel
  (\code{"g"}), Epanechnikov (\code{"e"}), uniform (\code{"u"}),
  quartic (\code{"q"}). The default (gaussian kernel) is strongly advised.}
  \item{ddlmini}{The number of eigenvalues (numerically) equals to 1.}
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
}

\value{
Returns the values of RMSE and MAP for each
value of the grid \code{K}. \code{Inf} are returned if the iteration leads
to a smoother with a df bigger than \code{ddlmaxi}.}

\references{
  Cornillon, P.-A.; Hengartner, N.; Jegou, N. and Matzner-Lober, E. (2012)
  Iterative bias reduction: a comparative study.
\emph{Statistics and Computing}. Doi: 10.1007/s11222-012-9346-4
}
\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\seealso{\code{\link{ibr}}}


\keyword{smooth}
\keyword{multivariate}
