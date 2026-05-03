predict.npregress <- function(object, newdata, interval= c("none", "confidence", "prediction"), deriv=FALSE, ...) {
  interval <- match.arg(interval)
  if ((interval == "prediction")|(interval == "confidence")) {
    warning("Interval for prediction/confidence is not implemented yet\n")
  }
    x <- object$call$x
    y <- object$call$y
  if (missing(newdata) || is.null(newdata)) {
    if (!deriv) return(object$fitted) else newdata <- x
  } else {
    if (any(is.na(newdata))) stop("NA's in newdata\n")
    if (!is.numeric(newdata)&(is.data.frame(newdata))) {
      newdata <- newdata[,1]
      if (!is.numeric(newdata)) stop("first column of data-frame is not numeric\n")
    }
    if (is.matrix(newdata)) {
      newdata <- as.vector(newdata)
    }
    if (!is.numeric(newdata)) stop("newdata must be a numeric vector (or a data-frame with first column of numeric type)\n")
  }
  kern <- c("g", "e", "q", "u")
  kernelint <- which(object$call$kernel==kern)
  ## autre methode
  if (object$call$degree==0) {
    methode <- "reg"
    prov <- .Call(ibr_npreg,as.double(x),as.double(y),as.double(newdata),as.double(object$bandwidth),as.integer(kernelint))
    deriv <- FALSE
  }
  if (object$call$degree==1) {
    prov <- .Call(ibr_npregpol,as.double(x),as.double(y),as.double(newdata),as.double(object$bandwidth),as.integer(kernelint))
  }
  if (!deriv) {
    Yres <- prov[[1]]
  } else {
    Yres <- list(yhat=prov[[1]],deriv=prov[[2]])  
  }
  if (object$call$degree>1) stop("Not implemented. Please consider using KernSmooth or another library for degree greater or equal to 2\n")
  return(Yres)
}

