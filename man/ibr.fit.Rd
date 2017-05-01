\name{ibr.fit}

\alias{ibr.fit}

\title{Iterative bias reduction smoothing}

\description{
  Performs iterative bias reduction using kernel, thin plate
  splines, Duchon splines or low rank splines. 
Missing values are not allowed. This function is not intended to be used directly.
}

\usage{
ibr.fit(x, y, criterion="gcv", df=1.5, Kmin=1, Kmax=1e+06, smoother="k",
 kernel="g", rank=NULL, control.par=list(), cv.options=list())}

\arguments{
  \item{x}{A numeric matrix of explanatory variables, with \emph{n} rows
    and \emph{p} columns.}
  \item{y}{A numeric vector of variable to be explained of length \emph{n}.}
  \item{criterion}{A vector of string. If the number of iterations
    (\code{iter}) is missing or 
    \code{NULL} the number of iterations is chosen using the either one
  criterion (the first
  coordinate of \code{criterion}) or several (see component
  \code{criterion} of argument list \code{control.par}). The criteria available are GCV (default, \code{"gcv"}),
    AIC (\code{"aic"}), corrected AIC  (\code{"aicc"}),   BIC
    (\code{"bic"}), gMDL  (\code{"gmdl"}), map (\code{"map"}) or rmse
  (\code{"rmse"}). The last two are designed for cross-validation.}
\item{df}{A numeric vector of either length 1 or length equal to the
  number of columns of \code{x}. If \code{smoother="k"}, it indicates
  the  desired effective degree of
  freedom (trace) of the smoothing   matrix for
  each variable or for the initial smoother (see \code{contr.sp$dftotal}); \code{df} is repeated when the length of vector
  \code{df} is 1. If \code{smoother="tps"} or  \code{smoother="ds"}, the
  minimum df of splines is multiplied by \code{df}. This argument is useless if
  \code{bandwidth} is supplied (non null).}
\item{Kmin}{The minimum number of bias correction  iterations of the
  search grid considered by
  the model selection procedure for selecting the optimal number of iterations.}
\item{Kmax}{The maximum number of bias correction  iterations of the
  search grid considered by
  the model selection procedure for selecting the optimal number of iterations.}
\item{smoother}{Character string which allows to choose between thin plate
  splines \code{"tps"}, Duchon
  splines \code{"tps"} (see Duchon, 1977) or kernel (\code{"k"}). }
\item{kernel}{Character string which allows to choose between gaussian kernel
  (\code{"g"}), Epanechnikov (\code{"e"}), uniform (\code{"u"}),
  quartic (\code{"q"}). The default (gaussian kernel) is strongly advised.}
\item{rank}{Numeric value that control the rank of low rank splines
  (denoted as \code{k} in mgcv package ; see also \link[mgcv]{choose.k}
  for further details or \link{gam} for another smoothing approach with
  reduced rank smoother.}

\item{control.par}{A named list that control optional parameters. The
  components are \code{bandwidth} (default to NULL), \code{iter}
  (default to NULL), \code{really.big} (default to \code{FALSE}),
  \code{dftobwitmax} (default to 1000), \code{exhaustive} (default to
  \code{FALSE}),\code{m} (default to NULL), ,\code{s} (default to NULL),
  \code{dftotal} (default to \code{FALSE}), \code{accuracy} (default to
  0.01), \code{ddlmaxi} (default to 2n/3), \code{fraction} (default
  to \code{c(100, 200, 500, 1000, 5000, 10^4, 5e+04, 1e+05, 5e+05,
    1e+06)}), \code{scale} (default to \code{FALSE}),
  \code{criterion} (default to \code{"strict"}) and
  \code{aggregfun} (default to 10^(floor(log10(x[2]))+2)).

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
  freedom and for \code{smoother="tps"} or \code{smoother="ds"} that lambda is chosen such that
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
  optimal number of iteration on the grid \code{Kmin:Kmax} is
  performed. All criteria for all iterations in the same class (class
  one: GCV, AIC, corrected AIC, BIC, gMDL ; class two : MAP, RMSE) are
  returned in argument \code{allcrit}. If \code{FALSE} the minimum of
  criterion is searched using \code{\link{optimize}} between \code{Kmin}
  and \code{Kmax}.
  
  \code{m}: The order of derivatives for the penalty (for thin plate
  splines it is the order). This integer \emph{m} must verify
  2\emph{m}+2\emph{s}/\emph{d}>1, where \emph{d} is the number of
  explanatory variables. The default (for \code{smoother="tps"}) is to
  choose the order \emph{m} as the first integer such that
  2\emph{m}/\emph{d}>1, where \emph{d} is the number of explanatory
  variables. The default (for \code{smoother="ds"}) is to choose
  \emph{m}=2 (p
  seudo cubic splines).

    \code{s}: the power of weighting function. For thin plate splines
  \emph{s} is equal to 0. This real must be strictly smaller than \emph{d}/2
  (where \emph{d} is the number of explanatory  variables) and must
  verify 2\emph{m}+2\emph{s}/\emph{d}. To get pseudo-cubic splines (the default),
  choose \emph{m}=2 and \emph{s}=(\emph{d}-1)/2 (See Duchon, 1977).the order of thin plate splines. This integer \emph{m} must verifies
  2\emph{m}/\emph{d}>1, where \emph{d} is the number of explanatory
  variables. 

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

  \code{scale}: boolean. If \code{TRUE} \code{x} is scaled (using
  \code{\link{scale}}); default to \code{FALSE}.

  \code{criterion} Character string. Possible choices are \code{strict},
  \code{aggregation} or \code{recalc}. \code{strict}
  allows to select the number of iterations according to
  the first coordinate of argument \code{criterion}.
  \code{aggregation}
  allows to select the number of iterations by applying the
  function \code{control.par$aggregfun} to the number of iterations
  selected by all the criteria chosen in argument \code{criterion}.
    \code{recalc}
    allows to select the number of iterations by first calculating the
    optimal number of the second coordinate of argument
  \code{criterion}, then applying the function
  \code{control.par$aggregfun} (to add some number to
  it) resulting in a new \code{Kmax} and then doing the optimal selction
  between  \code{Kmin} and this new \code{Kmax} using the first coordinate of argument
  \code{criterion}.
  ; default to \code{strict}.

  \code{aggregfun} function to be applied when
  \code{control.par$criterion} is either \code{recalc} or
  \code{aggregation}. 
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
Returns a list including:
  \item{beta}{Vector of coefficients.}
  \item{residuals}{Vector of residuals.}
  \item{fitted}{Vector of fitted values.}
  \item{iter}{The number of iterations used.}
  \item{initialdf}{The initial effective degree of freedom of the pilot (or base) smoother.}
  \item{finaldf}{The effective degree of freedom of the iterated bias reduction
    smoother at the \code{iter} iterations.}
  \item{bandwidth}{Vector of bandwith for each explanatory variable}
  \item{call}{The matched call}
  \item{parcall}{A list containing several components: 
    \code{p} contains the number of explanatory variables and \code{m}
    the order of the splines (if relevant), \code{s}
    the power of weights, \code{scaled} boolean which is \code{TRUE}
    when explanatory variables are scaled, \code{mean} mean of explanatory
    variables if \code{scaled=TRUE}, \code{sd} standard deviation of
    explanatory variables if \code{scaled=TRUE},  \code{critmethod} that indicates the method chosen
    for criteria \code{strict},
    \code{rank} the rank of low rank splines if relevant,
    \code{criterion} the chosen criterion,
    \code{smoother} the chosen smoother,
    \code{kernel} the chosen kernel,
    \code{smoothobject} the smoothobject returned by
    \link{smoothCon},
    \code{exhaustive} a boolean which indicates if an exhaustive
    search was chosen}
  \item{criteria}{Value
    of the chosen criterion at the given iteration, \code{NA} is
    returned when aggregation of criteria is chosen (see component
    \code{criterion} of list \code{control.par}). If the number of iterations
    \code{iter} is given by the user, \code{NULL} is returned}
  \item{alliter}{Numeric vector giving all the optimal number of iterations
    selected by  the chosen criteria.}
  \item{allcriteria}{either a list containing all the criteria evaluated on the
    grid \code{Kmin:Kmax} (along with the effective degree of freedom of the
    smoother and the sigma squared on this grid) if an exhaustive search is chosen (see the
    value of function 
    \code{\link{iterchoiceAe}} or \code{\link{iterchoiceS1e}}) 
     or all the values
    of criteria at the given optimal iteration if a non exhaustive
    search is chosen (see also \code{exhaustive} component of list
    \code{control.par}).}
}

\author{Pierre-Andre Cornillon, Nicolas Hengartner  and Eric Matzner-Lober.}

\references{
  Cornillon, P.-A.; Hengartner, N.; Jegou, N. and Matzner-Lober, E. (2012)
  Iterative bias reduction: a comparative study.
\emph{Statistics and Computing}, \emph{23}, 777-791.

  Cornillon, P.-A.; Hengartner, N. and Matzner-Lober, E. (2013)
  Recursive bias estimation for multivariate regression smoothers Recursive 
bias estimation for multivariate regression smoothers.
\emph{ESAIM: Probability and Statistics}, \emph{18}, 483-502.

  Cornillon, P.-A.; Hengartner, N. and Matzner-Lober, E. (2017)
  Iterative Bias Reduction Multivariate Smoothing in R: The ibr Package.
\emph{Journal of Statistical Software}, \emph{77}, 1--26.

  Wood, S.N. (2003) Thin plate regression
  splines. \emph{J. R. Statist. Soc. B}, \emph{65}, 95-114.
}

\seealso{\code{\link{ibr}}, \code{\link{predict.ibr}}, \code{\link{summary.ibr}}, \link[mgcv]{gam}}


\keyword{smooth}
\keyword{multivariate}
