\name{lrsmoother}

\alias{lrsmoother}

\title{Evaluate the lowrank spline}

\description{The function evaluates all the features needed for a
  lowrank spline smoothing. This function is not intended to be used directly.}

\usage{
lrsmoother(x,bs,listvarx,lambda,m,s,rank)
}

\arguments{
  \item{x}{Matrix of explanatory variables, size n,p.}
 \item{bs}{The type rank of lowrank splines: \code{tps} or \code{ds}.}
 \item{listvarx}{The vector of the names of explanatory variables}
  \item{lambda}{The smoothness coefficient lambda for thin plate splines of
  order \code{m}.}
  \item{m}{The order of derivatives for the penalty (for thin plate
  splines it is the order). This integer \emph{m} must verify
  2\emph{m}+2\emph{s}/\emph{d}>1, where \emph{d} is the number of explanatory
  variables.} 
  \item{s}{The power of weighting function. For thin plate splines
  \emph{s} is equal to 0. This real must be strictly smaller than \emph{d}/2
  (where \emph{d} is the number of explanatory  variables) and must
  verify 2\emph{m}+2\emph{s}/\emph{d}. To get pseudo-cubic splines,
  choose \emph{m}=2 and \emph{s}=(\emph{d}-1)/2 (See Duchon, 1977).} 
 \item{rank}{The rank of lowrank splines.}
 }
\details{see the reference for detailed explanation of the matrix matrix
  R^{-1}U (see reference) and \link{smoothCon} for the definition of smoothobject}
 
\references{
Duchon, J. (1977) Splines minimizing rotation-invariant semi-norms in
Solobev spaces. in W. Shemp and K. Zeller (eds) \emph{Construction theory of
functions of several variables}, 85-100, Springer, Berlin. 

Wood, S.N. (2003) Thin plate regression
  splines. \emph{J. R. Statist. Soc. B}, \emph{65}, 95-114.
}
 
\value{ Returns a list containing the smoothing matrix eigenvectors and eigenvalues
  \code{vectors} and \code{values},  and one
matrix denoted \code{Rm1U} and one smoothobject \code{smoothobject}.}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober}

\seealso{\code{\link{ibr}}}


\keyword{smooth}
