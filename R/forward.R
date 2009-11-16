forward <- function(X,Y,criterion="gcv",df=1.5,Kmin=1,Kmax=10000,smoother="k",kernel="g",control.par=list(),cv.options=list()) {
  garde <- NULL
  ptot <- ncol(X)
  possible <- 1:ptot
  resultat <- matrix(Inf,ptot,ptot)
  letop <- Inf
  bandwidth <- NULL
  if (!((any(names(control.par)=="bandwidth")&&is.numeric(control.par$bandwidth))|((any(names(control.par)=="dftotal"))&&control.par$dftotal))&(smoother=="k")) {
    if (any(names(control.par)=="dftobwitmax")) dftobwitmax <- control.par$dftobwitmax else dftobwitmax <- 1000
    if (!is.numeric(dftobwitmax)) stop("argument dftobwitmax of list control.par is not numeric\n")
    if (dftobwitmax<1) stop("argument dftobwitmax of list control.par is less than 1\n")
    bandwidth <- bwchoice(X,df,kernel,dftobwitmax)
  }
  for (iter in 1:ptot) {
      for (i in possible) {
        garden <- sort(c(garde,i))
        if (is.null(bandwidth)) {
          res3 <- ibr(X[,garden,drop=FALSE],Y,criterion,df,Kmin,Kmax,smoother,kernel,control.par,cv.options)
        } else {
          res3 <- ibr(X[,garden,drop=FALSE],Y,criterion,df,Kmin,Kmax,smoother,kernel,control.par=c(control.par,list(bandwidth=bandwidth[garden])),cv.options)
        }
        if (any(names(control.par)=="exhaustive")&&control.par$exhaustive) {
          allcrit <- switch(criterion,rmse=res3$criteria$rmse,map=res3$criteria$map,aic=res3$criteria$aic,aicc=res3$criteria$aicc,gcv=res3$criteria$gcv,bic=res3$criteria$bic,gmdl=res3$criteria$gmdl)
          allcrit <- allcrit[res3$iter]
        }
        else allcrit <- res3$criteria
        resultat[iter,i] <- allcrit
      }
      mini <- which.min(resultat[iter,])
      letopnew <-min(resultat[iter,])
      if (letopnew>letop) {
        resultat <- resultat[1:(iter-1),]
        if (iter>1) {
          class(resultat) <- c("forwardibr", "matrix")
        } else {
            class(resultat) <- c("forwardibr", "vector")
          }
        return(resultat) } else letop <- letopnew
      coordonnee <- which(possible==mini)
      garde <- sort(c(garde,possible[coordonnee]))
      possible <- possible[-coordonnee]
    }
  if (iter>1) {
    class(resultat) <- c("forwardibr", "matrix")
  } else {
    class(resultat) <- c("forwardibr", "vector")
  }
  return(resultat)
}

