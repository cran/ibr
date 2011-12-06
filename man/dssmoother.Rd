\name{dssmoother}

\alias{dssmoother}

\title{Evaluate the smoothing matrix, the radial basis matrix, the
  polynomial matrix and their associated coefficients}

\description{The function evaluates the smoothing matrix \code{H}, the
  matrices \emph{Q} and \emph{S} and their associated
  coefficients \code{c} and \code{s}. This function is not intended to be used directly.}

\usage{
dssmoother(X,Y=NULL,lambda,m,s)
}

\arguments{
  \item{X}{Matrix of explanatory variables, size n,p.}
  \item{Y}{Vector of  response variable. If null, only the smoothing
  matrix is returned.}
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
 }
\details{see the reference for detailed explanation of \emph{Q} (the
  semi kernel or radial basis) and
  \emph{S} (the polynomial null space). Both matrices are computed using
  \pkg{fields}.}
 
\references{
Duchon, J. (1977) Splines minimizing rotation-invariant semi-norms in
Solobev spaces. in W. Shemp and K. Zeller (eds) \emph{Construction theory of
functions of several variables}, 85-100, Springer, Berlin. 

C. Gu (2002) \emph{Smoothing spline anova models}. New York:
Springer-Verlag.
}
 
\value{ Returns a list containing the smoothing matrix \code{H}, and two
matrices denoted \code{Sgu} (for null space) and \code{Qgu}.}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober}

\seealso{\code{\link{ibr}}, \code{\link[fields]{fields.mkpoly}}, \code{\link[fields]{Rad.cov}}}


\keyword{smooth}
