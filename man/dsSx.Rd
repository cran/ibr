\name{dsSx}

\alias{dsSx}

\title{Evaluate the smoothing matrix at any point}

\description{The function evaluates the matrix \emph{Q} and \emph{S} related to
  the explanatory variables \eqn{X} at any points. This function is not intended to be used directly.
}

\usage{
dsSx(X,Xetoile,m=2,s=0)
}

\arguments{
  \item{X}{Matrix of explanatory variables, size n,p.}
  \item{Xetoile}{Matrix of new observations with the same number of
    variables as \eqn{X}, size m,p.}
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
  semi kernel) and
  \emph{S} (the polynomial null space).}
 
\references{
Duchon, J. (1977) Splines minimizing rotation-invariant semi-norms in
Solobev spaces. in W. Shemp and K. Zeller (eds) \emph{Construction theory of
  functions of several variables}, 85-100, Springer, Berlin.

C. Gu (2002) \emph{Smoothing spline anova models}. New York:
Springer-Verlag.
}
 
\value{
Returns a list containing two matrices denoted \code{Sgu} (for null
space) and \code{Qgu}}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober}

\seealso{\code{\link{ibr}}}


\keyword{smooth}
