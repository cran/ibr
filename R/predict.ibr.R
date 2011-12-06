predict.ibr <- function(object, newdata, interval= c("none", "confidence", "prediction"), ...) {
  interval <- match.arg(interval)
    if ((interval == "prediction")|(interval == "confidence")) {
      warning("Interval for prediction/confidence is not implemented yet\n")
    }
  if (missing(newdata) || is.null(newdata)) {
    Yres <- object$fitted
  } else {
    if (any(is.na(newdata))) stop("NA's in newdata\n")
    if (!is.matrix(newdata)) {
      newdata <- data.matrix(newdata)
      warning("newdata is coerced to matrix by data.matrix function\n")
    }
    if (object$call$scaled) {
      newdata  <- scale(newdata,center=object$call$mean,scale=object$call$sd)
    }
    if (ncol(newdata)!=object$call$p) stop("the number of variables in new dataset is different from the one in the original one\n")
    x <- object$call$x
    if (object$call$smoother=="k") {
      SSx <- kernelSx(kernelx=object$call$kernel,X=x,bx=object$bandwidth,newdata)
      Yres <- as.vector(SSx%*%object$beta)
    }
    if ((object$call$smoother=="ds")|(object$call$smoother=="tps")) {
      SSx <- dsSx(x,newdata,m=object$call$m,s=object$call$s)
      Yres <- as.vector(SSx$Sgu%*%object$beta$d+SSx$Qgu%*%object$beta$c)
    }    
  }
  return(Yres)
}

