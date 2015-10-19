\name{betaS1lr}

\alias{betaS1lr}

\title{Coefficients for iterative bias reduction method. }

\description{The function evaluates the smoothing matrix \code{H}, the
  matrices \emph{Q} and \emph{S} and their associated
  coefficients \code{c} and \code{s}. This function is not intended to be used directly.}

\usage{
betaS1lr(n,U,tUy,eigenvaluesS1,ddlmini,k,lambda,rank,Rm1U,index0)
}

\arguments{
  \item{n}{The number of observations.}
  \item{U}{The the matrix of eigen vectors of the
  symmetric smoothing matrix \emph{S}.}
  \item{tUy}{The transpose of the matrix of eigen vectors of the
  symmetric smoothing matrix \emph{S} times the vector of observation \emph{y}.}
  \item{eigenvaluesS1}{Vector of the eigenvalues of the
  symmetric smoothing matrix \emph{S}.}
  \item{ddlmini}{The number of eigen values of \emph{S} equal to 1.}
  \item{k}{A numeric vector which give the number of iterations.}
 \item{lambda}{The smoothness coefficient lambda for thin plate splines of
  order \code{m}.}
 \item{rank}{The rank of lowrank splines.}
 \item{Rm1U}{matrix R^{-1}U (see reference).}
  \item{index0}{The index of the first eigen values of \emph{S}
  numerically equal to 0.}
 }
\details{See the reference for detailed explanation of \emph{Q} (the
  semi kernel or radial basis) and
  \emph{S} (the polynomial null space).}
 
\references{
Wood, S.N. (2003) Thin plate regression
  splines. \emph{J. R. Statist. Soc. B}, \emph{65}, 95-114.
} 
\value{
Returns \code{beta}}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober}

\seealso{\code{\link{ibr}}}


\keyword{smooth}
