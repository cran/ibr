\name{ibr}

\alias{ibr}
\alias{print.ibr}
\alias{residuals.ibr}

\title{Iterative bias reduction smoothing}

\description{
  Performs iterative bias reduction using kernel or thin plate
  splines. In the latter case, the order \emph{m} is chosen as the first integer
  such that 2\emph{m}/\emph{d}>1, where \emph{d} is the number of explanatory variables.
  \cr
Missing values are not allowed.
}

\usage{
ibr(x, y, criterion="gcv", df=1.5, Kmin=1, Kmax=10000, smoother="k",
 kernel="g", control.par=list(), cv.options=list())}

\arguments{
  \item{x}{A numeric matrix of explanatory variables, with \emph{n} rows
    and \emph{p} columns.}
  \item{y}{A numeric vector of variable to be explained of length \emph{n}.}
  \item{criterion}{Character string. If the number of iterations
    (\code{iter}) is missing or 
    \code{NULL} the number of iterations is chosen using
    \code{criterion}. The criteria available are GCV (default, \code{"gcv"}),
    AIC (\code{"aic"}), corrected AIC  (\code{"aicc"}),   BIC
    (\code{"bic"}), gMDL  (\code{"gmdl"}), map (\code{"map"}) or rmse
  (\code{"rmse"}). The last two are designed for cross-validation.}
\item{df}{A numeric vector of either length 1 or length equal to the
  number of columns of \code{x}. If \code{smoother="k"}, it indicates
  the  desired effective degree of
  freedom (trace) of the smoothing   matrix for
  each variable or for the initial smoother (see \code{contr.sp$dftotal}); \code{df} is repeated when the length of vector
  \code{df} is 1. If \code{smoother="tps"}, the minimum df of thin
  plate splines is multiplied by \code{df}. This argument is useless if
  \code{bandwidth} is supplied (non null).}
\item{Kmin}{The minimum number of bias correction  iterations of the
  search grid considered by
  the model selection procedure for selecting the optimal number of iterations.}
\item{Kmax}{The maximum number of bias correction  iterations of the
  search grid considered by
  the model selection procedure for selecting the optimal number of iterations.}
\item{smoother}{Character string which allows to choose between thin plate
  splines \code{"tps"} or kernel (\code{"k"}). }
\item{kernel}{Character string which allows to choose between gaussian kernel
  (\code{"g"}), Epanechnikov (\code{"e"}), uniform (\code{"u"}),
  quartic (\code{"q"}). The default (gaussian kernel) is strongly advised.}

\item{control.par}{A named list that control optional parameters. The
  components are \code{bandwidth} (default to NULL), \code{iter}
  (default to NULL), \code{really.big} (default to \code{FALSE}),
  \code{dftobwitmax} (default to 1000), \code{exhaustive} (default to
  \code{FALSE}),\code{m} (default to NULL), \code{dftotal} (default to
  \code{FALSE}), \code{accuracy} (default to 0.01), \code{ddlmaxi}
  (default to 2n/3) and \code{fraction} (default to \code{c(100, 200,
  500, 1000, 5000, 10^4, 5e+04, 1e+05, 5e+05, 1e+06)}).

  \code{bandwidth}: a vector of either length 1 or length equal to the
    number of columns of \code{x}. If \code{smoother="k"},
    it indicates the bandwidth used for
  each variable, bandwidth is repeated when the length of vector
  \code{bandwidth} is 1. If \code{smoother="tps"}, it indicates the
  amount of penalty (coefficient lambda).
  The default (missing) indicates, for \code{smoother="k"}, that
  bandwidth for each variable is
  chosen such that each univariate kernel
  smoother (for each explanatory variable) has \code{df} effective degrees of
  freedom and for \code{smoother="tps"} that lambda is chosen such that
  the df of the smoothing matrix is \code{df} times the minimum df.

  \code{iter}: the number of iterations. If null or missing, an optimal number of
  iterations is chosen from 
  the search grid (integer from \code{Kmin} to \code{Kmax}) to minimize the \code{criterion}.
 
  \code{really.big}: a boolean: if \code{TRUE} it overides the limitation
  at 500 observations. Expect long computation times if \code{TRUE}.
  
\code{dftobwitmax}: When bandwidth is chosen by specifying the effective
  degree
  of freedom (see \code{df}) a search is done by
  \code{\link{uniroot}}. This argument specifies the maximum number of iterations transmitted to \code{\link{uniroot}} function.
  
   \code{exhaustive}: boolean, if \code{TRUE} an exhaustive search of
  optimal number of iteration on the
  grid \code{Kmin:Kmax} is performed. If \code{FALSE} the minimum  of
  criterion is searched using \code{\link{optimize}} between \code{Kmin}
  and \code{Kmax}.
  
  \code{m}: the order of thin plate splines. This integer \emph{m} must verifies
  2\emph{m}/\emph{d}>1, where \emph{d} is the number of explanatory
  variables.  The missing default to choose the order \emph{m} as the first integer
  such that 2\emph{m}/\emph{d}>1, where \emph{d} is the number of
  explanatory variables (same for \code{NULL}). 
  
  \code{dftotal}: a boolean wich indicates when \code{FAlSE} that the
  argument \code{df} is the objective df for each univariate kernel (the
  default) calculated for each explanatory variable or for the overall
  (product) kernel, that is the base smoother (when \code{TRUE}).
  
  \code{accuracy}: tolerance when searching bandwidths which lead to a
  chosen overall intial df.
  
  \code{dfmaxi}: the maximum effective degree of freedom allowed for iterated
  biased reduction smoother. 
  
  \code{fraction}: the subdivision of interval \code{Kmin},\code{Kmax}
  if non exhaustive search is performed (see also \code{\link{iterchoiceA}} or \code{\link{iterchoiceS1}}). 
}

   \item{cv.options}{A named list which controls the way to do cross
  validation with component \code{bwchange},
  \code{ntest}, \code{ntrain}, \code{Kfold}, \code{type},
  \code{seed}, \code{method} and \code{npermut}. \code{bwchange} is a boolean (default to \code{FALSE})
  which indicates if bandwidth have to be recomputed each
  time. \code{ntest} is the number of observations in test set and
  \code{ntrain} is the number of observations in training set. Actually,
  only one of these is needed the other can be \code{NULL} or missing. \code{Kfold} a boolean or an integer. If
  \code{Kfold} is \code{TRUE} then the number of fold is deduced from
  \code{ntest} (or \code{ntrain}).  \code{type} is a character string in
  \code{random},\code{timeseries},\code{consecutive}, \code{interleaved}
  and give the type of segments.  \code{seed} controls the seed of
  random generator. \code{method} is either \code{"inmemory"} or
  \code{"outmemory"}; \code{"inmemory"} induces some calculations outside
  the loop saving computational time but leading to an increase of the required
  memory. \code{npermut} is the number of random draws.   If
  \code{cv.options} is \code{list()}, then component \code{ntest} is set to
  \code{floor(nrow(x)/10)}, \code{type} is random, \code{npermut} is 20
  and \code{method} is \code{"inmemory"}, and the other components are \code{NULL}}
}

\value{
Returns an object of class \code{ibr} which is a list including:
  \item{beta}{Vector of coefficients.}
  \item{residuals}{Vector of residuals.}
  \item{iter}{The number of iterations used.}
  \item{initialdf}{The initial effective degree of freedom of the pilot (or base) smoother.}
  \item{finaldf}{The effective degree of freedom of the iterated bias reduction
    smoother at the \code{iter} iterations.}
  \item{bandwidth}{Vector of bandwith for each explanatory variable}
  \item{call}{A list containing four components: \code{x} contains the
    initial explanatory variables, \code{y} contains the
    initial dependant variables, 
    \code{criterion} contains the chosen criterion, \code{kernel} the
    kernel, \code{p} the number of explanatory variables and \code{m}
    the order of the splines (if relevant).}
  \item{criteria}{either a list containing all the criteria evaluated on the
    grid \code{Kmin:Kmax} (along with the effective degree of freedom of the
    smoother and the sigma squared on this grid) if an exhaustive search is chosen (see the
    value of function 
    \code{\link{iterchoiceAe}} or \code{\link{iterchoiceS1e}}) or the value
    of the chosen criterion at the given iteration if a non exhaustive
    search is chosen (see \code{exhaustive}). If the number of iterations
    \code{iter} is given by the user \code{NULL} is returned}
}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\references{
Cornillon, P. A., Hengartner, N. and Matzner-Lober,  E. (2009) Recursive
Bias Estimation for high dimensional regression smoothers. \emph{submitted}.}

\seealso{\code{\link{predict.ibr}}, \code{\link{summary.ibr}}}

\examples{
f <- function(x, y) { .75*exp(-((9*x-2)^2 + (9*y-2)^2)/4) +
                      .75*exp(-((9*x+1)^2/49 + (9*y+1)^2/10)) +
                      .50*exp(-((9*x-7)^2 + (9*y-3)^2)/4) -
                      .20*exp(-((9*x-4)^2 + (9*y-7)^2)) }
# define a (fine) x-y grid and calculate the function values on the grid
ngrid <- 50; xf <- seq(0,1, length=ngrid+2)[-c(1,ngrid+2)]
yf <- xf ; zf <- outer(xf, yf, f)
grid <- cbind(rep(xf, ngrid), rep(xf, rep(ngrid, ngrid)))
persp(xf, yf, zf, theta=130, phi=20, expand=0.45,main="True Function")
#generate a data set with function f and noise to signal ratio 5
noise <- .2 ; N <- 100 
xr <- seq(0.05,0.95,by=0.1) ; yr <- xr ; zr <- outer(xr,yr,f) ; set.seed(25)
std <- sqrt(noise*var(as.vector(zr))) ; noise <- rnorm(length(zr),0,std)
Z <- zr + matrix(noise,sqrt(N),sqrt(N))
# transpose the data to a column format 
xc <- rep(xr, sqrt(N)) ; yc <- rep(yr, rep(sqrt(N),sqrt(N)))
X <- cbind(xc, yc) ; Zc <- as.vector(Z)
# fit by thin plate splines (of order 2) ibr
res.ibr <- ibr(X,Zc,df=1.1,smoother="tps")
fit <- matrix(predict(res.ibr,grid),ngrid,ngrid)
persp(xf, yf, fit ,theta=130,phi=20,expand=0.45,main="Fit",zlab="fit")

\dontrun{
data(ozone, package = "ibr")
res.ibr <- ibr(ozone[,-1],ozone[,1],df=1.1)
summary(res.ibr)
predict(res.ibr)}
}

\keyword{smooth}
\keyword{multivariate}
