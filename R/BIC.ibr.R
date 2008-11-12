BIC <- function(object, ...) UseMethod("BIC")
AICc <- function(object, ...) UseMethod("AICc")
GCV <- function(object, ...) UseMethod("GCV")

BIC.ibr <- function(object, ...) {
  r <- object$residuals
  n <- length(r)
  stderr <- sqrt(sum(r^2)/(n-object$finaldf))
  return(log(stderr^2)+log(n)*object$finaldf/n)
}
GCV.ibr <- function(object, ...) {
  r <- object$residuals
  n <- length(r)
  stderr <- sqrt(sum(r^2)/(n-object$finaldf))
  return(log(stderr^2)-2*log(1-object$finaldf/n))
}
AICc.ibr <- function(object, ...) {
  r <- object$residuals
  n <- length(r)
  stderr <- sqrt(sum(r^2)/(n-object$finaldf))
  return(log(stderr^2)+1+(2*(object$finaldf+1))/(n-object$finaldf-2))
}
