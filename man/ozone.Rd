\name{ozone}
\docType{data}
\alias{ozone}
\title{Los Angeles ozone pollution data, 1976.}

\description{
Los Angeles ozone pollution data, 1976.  
We deleted from the original data, the first 3 columns which were the  
\code{Month}, \code{Day of the month} and \code{Day of the week}. Each
observation is one day, so there is 366 rows.
The \code{ozone} data is a matrix with 9 columns.
}

\format{
This data set is a matrix containing the following columns: 
\tabular{rlll}{
         [,1] \tab Ozone   \tab numeric \tab Daily maximum
	 one-hour-average ozone reading 
	 (parts per million) at Upland, CA.\cr
         [,2] \tab Pressure.Vand \tab numeric \tab 500 millibar pressure height (m) measured at Vandenberg AFB.\cr
         [,3] \tab Wind    \tab numeric \tab Wind speed (mph) at Los Angeles International Airport (LAX).\cr
         [,4] \tab Humidity \tab numeric \tab Humidity in percentage at LAX.\cr
         [,5] \tab Temp.Sand \tab numeric \tab Temperature (degrees F) measured at Sandburg, CA.\cr
         [,6] \tab Inv.Base.height \tab numeric \tab Inversion base height (feet) at LAX.\cr
         [,7] \tab Pressure.Grad   \tab numeric \tab Pressure gradient (mm Hg) from LAX to Daggett, CA.\cr
         [,8] \tab Inv.Base.Temp   \tab numeric \tab Inversion base temperature (degrees F) at LAX.\cr
         [,9] \tab Visilibity     \tab numeric \tab Visibility (miles) measured at LAX.
       }
     }
\source{
Leo Breiman, Department of Statistics, UC Berkeley.  Data used in
     Breiman, L. and  Friedman, J. H. (1985). Estimating optimal
     transformations for multiple regression and correlation,
     \emph{Journal of American Statistical Association}, \bold{80},
     580--598.
}
\seealso{ \code{\link{ibr}}}
\keyword{datasets}
