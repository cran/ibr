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
  ## autre methode
  if (object$call$degree==0) {
    methode <- "reg"
    if (object$call$kernel=="g") {
    prov <- .C("regg",as.double(x),as.integer(length(x)),as.double(y),as.double(object$bandwidth),as.double(newdata),as.integer(length(newdata)),double(length(newdata)),double(1))
    }
    if (object$call$kernel=="q") {
    prov <- .C("regq",as.double(x),as.integer(length(x)),as.double(y),as.double(object$bandwidth),as.double(newdata),as.integer(length(newdata)),double(length(newdata)),double(1))
    }
    if (object$call$kernel=="e") {
    prov <- .C("rege",as.double(x),as.integer(length(x)),as.double(y),as.double(object$bandwidth),as.double(newdata),as.integer(length(newdata)),double(length(newdata)),double(1))
    }
    if (object$call$kernel=="u") {
    prov <- .C("regu",as.double(x),as.integer(length(x)),as.double(y),as.double(object$bandwidth),as.double(newdata),as.integer(length(newdata)),double(length(newdata)),double(1))
    }
    deriv <- FALSE
  }
  if (object$call$degree==1) {
    if (object$call$kernel=="g") {
    prov <- .C("regpolg",as.double(x),as.integer(length(x)),as.double(y),as.double(object$bandwidth),as.double(newdata),as.integer(length(newdata)),double(length(newdata)),double(1),double(length(newdata)))
    }
    if (object$call$kernel=="q") {
    prov <- .C("regpolq",as.double(x),as.integer(length(x)),as.double(y),as.double(object$bandwidth),as.double(newdata),as.integer(length(newdata)),double(length(newdata)),double(1),double(length(newdata)))
    }
    if (object$call$kernel=="e") {
    prov <- .C("regpole",as.double(x),as.integer(length(x)),as.double(y),as.double(object$bandwidth),as.double(newdata),as.integer(length(newdata)),double(length(newdata)),double(1),double(length(newdata)))
    }
    if (object$call$kernel=="u") {
    prov <- .C("regpolu",as.double(x),as.integer(length(x)),as.double(y),as.double(object$bandwidth),as.double(newdata),as.integer(length(newdata)),double(length(newdata)),double(1),double(length(newdata)))
    }
  }
  if (!deriv) {
    Yres <- prov[[7]]
  } else {
    Yres <- list(yhat=prov[[7]],deriv=prov[[9]])  
  }
  if (object$call$degree>1) stop("Not implemented. Please consider using KernSmooth or another library for degree greater or equal to 2\n")
  return(Yres)
}

