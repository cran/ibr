\name{cvobs}

\alias{cvobs}

\title{Selection of the number of iterations for iterative bias reduction smoothers}

\description{
The function \code{cvobs} gives the index of observations in each test set. This function is not intended to be used directly.
}

\usage{
cvobs(n,ntest,ntrain,Kfold,type=
c("random", "timeseries", "consecutive", "interleaved"), npermut, seed)}

\arguments{
  \item{n}{The total number of observations.}
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
}

\value{
Returns a list with in each component the index of observations to be
used as a test set.}

\references{
Cornillon, P. A., Hengartner, N. and Matzner-Lober,  E. (2009) Recursive
Bias Estimation for high dimensional regression smoothers. \emph{submitted}.}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\seealso{\code{\link{ibr}}}


\keyword{smooth}
\keyword{multivariate}
