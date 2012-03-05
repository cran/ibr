\name{betaS1}

\alias{betaS1}

\title{Coefficients for iterative bias reduction method. }

\description{The function evaluates the smoothing matrix \code{H}, the
  matrices \emph{Q} and \emph{S} and their associated
  coefficients \code{c} and \code{s}. This function is not intended to be used directly.}

\usage{
betaS1(n,U,tUy,eigenvaluesS1,ddlmini,k,lambda,Sgu,Qgu,index0)
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
 \item{Sgu}{The matrix of the polynomial null space \emph{S}.}
 \item{Qgu}{The matrix of the semi kernel (or radial basis) \emph{Q}.}
  \item{index0}{The index of the first eigen values of \emph{S}
  numerically equal to 0.}
 }
\details{See the reference for detailed explanation of \emph{Q} (the
  semi kernel or radial basis) and
  \emph{S} (the polynomial null space).}
 
\references{
C. Gu (2002) \emph{Smoothing spline anova models}. New York: Springer-Verlag.}
 
\value{
Returns a list containing of coefficients for the null space \code{dgub}
and the semi-kernel \code{cgub}}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober}

\seealso{\code{\link{ibr}}}


\keyword{smooth}
