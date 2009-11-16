\name{tpsSx}

\alias{tpsSx}

\title{Evaluate the smoothing matrix at any point}

\description{The function evaluates the matrix \emph{Q} and \emph{S} related to
  the explanatory variables \eqn{X} at any points. This function is not intended to be used directly.
}

\usage{
tpsSx(X,Xetoile,lambda,m=2)
}

\arguments{
  \item{X}{Matrix of explanatory variables, size n,p.}
  \item{Xetoile}{Matrix of new observations with the same number of
    variables as \eqn{X}, size m,p.}
 \item{lambda}{Coefficient lambda that control smoothness.}
 \item{m}{Order of thin plate splines.}
 }
\details{see the reference for detailed explanation of \emph{Q} (the
  semi kernel) and
  \emph{S} (the polynomial null space). Both matrices are computed using
  \pkg{fields}.}
 
\references{
C. Gu (2002) \emph{Smoothing spline anova models}. New York: Springer-Verlag.}
 
\value{
Returns a list containing two matrices denoted \code{Sgu} (for null
space) and \code{Qgu}}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober}

\seealso{\code{\link{ibr}}, \code{\link[fields]{fields.mkpoly}}, \code{\link[fields]{Rad.cov}}}


\keyword{smooth}
