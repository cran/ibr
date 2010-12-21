summary.ibr <- function(object, criteria="call", ...) {
  r <- object$residuals
  x <- object$call$x
  y <- object$call$y
  n <- length(r)
  stderr <- sqrt(sum(r^2)/(n-object$finaldf))
  sigma2 <- stderr^2
  if (any(criteria=="call")) {
    criteria <- object$call$criterion
    anscrit <- NULL
  } else {
    crit <-c("aic","aicc","gcv","bic","gmdl")
    if (all(!(criteria%in%crit))) stop(paste("criteria are:",crit,"\n"))
    criteria <- criteria[criteria%in%crit]
    anscrit <- NULL
  }
  if (any(criteria=="gcv"))  anscrit <- c(anscrit,log(sigma2)-2*log(1-object$finaldf/n))
  if (any(criteria=="aic"))  anscrit <- c(anscrit,log(sigma2)+2*object$finaldf/n)
  if (any(criteria=="aicc"))  anscrit <- c(anscrit,log(sigma2)+1+(2*(object$finaldf+1))/(n-object$finaldf-2))
  if (any(criteria=="bic"))  anscrit <- c(anscrit,log(sigma2) + log(n)*(object$finaldf)/n)
  if (any(criteria=="gmdl")) {
    Sbul <-   n*sigma2/(n-object$finaldf)
  anscrit <- c(anscrit,log(Sbul)+object$finaldf/n*log((sum(y^2)-n*sigma2)/(object$finaldf*Sbul)))
  }
  if ((criteria!="user")&(!is.null(anscrit))) {
  names(anscrit) <- criteria
} else {
  anscrit <- "No Informative Criterion"
  names(anscrit) <- criteria
}
  ans <- list(residuals=r,Std.Error=stderr,Initial.Df=object$initialdf,
              Final.Df=object$finaldf,Resid.Df=n-object$finaldf,criteria=anscrit,
              kernel=object$call$kernel, iter=object$iter,crit4iter=object$call$criterion,
              bandwidth=object$bandwidth, smoother=object$call$smoother,m=object$call$m)
  if (object$call$smoother=="k") {
    if (is.null(colnames(x))) {
      names(ans$bandwidth) <- paste("X",1:object$call$p,sep="")
    } else {
      names(ans$bandwidth) <- paste(colnames(x),1:object$call$p,sep="")
    }
  }
  class(ans) <- "summary.ibr"
  ans
}
