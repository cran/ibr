\name{tpssmoother}

\alias{tpssmoother}

\title{Evaluate the smoothing matrix, the radial basis matrix, the
  polynomial matrix and their associated coefficients}

\description{The function evaluates the smoothing matrix \code{H}, the
  matrices \emph{Q} and \emph{S} and their associated
  coefficients \code{c} and \code{s}. This function is not intended to be used directly.}

\usage{
tpssmoother(X,Y=NULL,lambda,m=2)
}

\arguments{
  \item{X}{Matrix of explanatory variables, size n,p.}
  \item{Y}{Vector of  response variable. If null, only the smoothing
  matrix is returned.}
  \item{lambda}{The smoothness coefficient lambda for thin plate splines of
  order \code{m}.}
 \item{m}{Order of thin plate splines.}
 }
\details{see the reference for detailed explanation of \emph{Q} (the
  semi kernel or radial basis) and
  \emph{S} (the polynomial null space). Both matrices are computed using
  \pkg{fields}.}
 
\references{
C. Gu (2002) \emph{Smoothing spline anova models}. New York: Springer-Verlag.}
 
\value{ Returns a list containing the smoothing matrix \code{H}, and two
matrices denoted \code{Sgu} (for null space) and \code{Qgu}.}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober}

\seealso{\code{\link{ibr}}, \code{\link[fields]{fields.mkpoly}}, \code{\link[fields]{Rad.cov}}}


\keyword{smooth}
