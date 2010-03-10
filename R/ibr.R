ibr <- function(x,y,criterion="gcv",df=1.5,Kmin=1,Kmax=10000,smoother="k",kernel="g",control.par=list(),cv.options=list()) {
  crit <-c("aic","aicc","gcv","bic","gmdl","rmse","map")
  criterion <- match.arg(criterion,crit)
  smoothertable <- c("k","tps")
  smoother <- match.arg(smoother,smoothertable)
  if (!is.matrix(x)) {
    x <- data.matrix(x)
    warning("x is coerced to matrix by data.matrix function\n")
  }
  n <- nrow(x)
  p <- ncol(x)
  if (length(y)!=n) stop("numbers of observations in x and y do not match\n")
  contr.sp <- list(bandwidth=NULL,iter=NULL,really.big=FALSE,dftobwitmax=1000,exhaustive=FALSE,m=NULL,dftotal=FALSE,accuracy=0.01,dfmaxi=2*n/3,fraction=c(100, 200, 500, 1000, 5000,10^4,5e+04,1e+05,5e+05,1e+06))
  contr.sp[(names(control.par))] <- control.par
  if (!(is.logical(contr.sp$dftotal))) stop("contr.sp$dftotal must be logical\n")
  if (is.null(contr.sp$m)) contr.sp$m <- floor(p/2)+1 else {
    if (contr.sp$m<=(p/2)) stop("order of thin plate splines is invalid\n")
  }
  if ((!is.null(contr.sp$bandwidth))&(!is.numeric(contr.sp$bandwidth) || (contr.sp$bandwidth<0))) stop("invalid bandwidth\n")
  if ((!is.null(contr.sp$iter))&(!is.numeric(contr.sp$iter) || (contr.sp$iter<0) || (floor(contr.sp$iter)!=contr.sp$iter))) stop("invalid number of iterations\n")
  iter <- contr.sp$iter
   if ((contr.sp$dfmaxi<=0)|(contr.sp$dfmaxi>n)) stop("invalid dfmaxi\n")
  if (smoother=="tps") {
    if (!is.numeric(contr.sp$m) || (contr.sp$m<0) || (floor(contr.sp$m)!=contr.sp$m)) stop("invalid spline order\n")
  } else contr.sp$m <- NULL
  if (criterion%in%c("rmse","map")) {
    cv <- list(bwchange=FALSE,ntest=floor(nrow(x)/10),ntrain=NULL,Kfold=FALSE,type="random",seed=NULL,npermut=20)
    cv[(names(cv.options))] <- cv.options
    if (!all(sapply(cv[1],is.logical))) stop("invalid cv$bwchange or cv$Kfold: must be logical\n")
    if (!all(sapply(cv[c(2,3,6,8)], FUN=function(x) is.numeric(x)||is.null(x)))) stop("invalid cv parameters: must be numeric or NULL\n")
    if (any(names(cv.options)=="ntrain")) cv$ntest <- NULL
  } else cv <- NULL
  if ((n>500)&(! contr.sp$really.big)) stop("number of observations is greater than 500, set control.par$really.big to TRUE to you do want to do calculations (but computation time -eigen decomposition- could be prohibitive)\n")
  if (smoother=="k") {
    m <- NULL
    if (!is.null(contr.sp$bandwidth)) {
      if (length(contr.sp$bandwidth)==1) {
        bandwidth <- rep(contr.sp$bandwidth,p)
      } else {
        if (length(contr.sp$bandwidth)!=p)  stop(paste("the length of bandwidth vector have to be",p,"or 1\n"))
        bandwidth <- contr.sp$bandwidth
      }
      listeA <- calcA(X = x, bx = bandwidth, kernelx = kernel)
    } else {
      if (kernel=="g") {
        if (df<=1) stop("degree of freedom should be greater than 1\n")   
        departbw <- apply(x,2,FUN=function(z) 3*abs(diff(range(z))))
        calculus <- .C("gaustotal",
                      as.double(departbw),
       if (contr.sp$dftotal) as.double(rep(1e-10,p)) else
                       as.double(1e-10),
                 as.double(x),as.integer(n),as.integer(p),
                 as.double(.Machine$double.eps^0.25),
                 maxit=as.integer(contr.sp$dftobwitmax),
                 as.double(df),as.integer(contr.sp$dftotal),
                 A=double(n^2),Ddemi=double(n),df=double(1),
                 bandwidth=double(p),PACKAGE="ibr",DUP=FALSE)
        if (calculus$maxit==-1) warning("failed to find appropriate bandwidth. Try to increase dftobwitmax argument of control.par list ?\n")
        listeA <- list(A=matrix(calculus$A,n,n),Ddemi=calculus$Ddemi,df=calculus$df)
        bandwidth <- calculus$bandwidth
      rm(calculus)
      } else {
        if (contr.sp$dftotal) {
          dfobjectif <- df
          dfmini <- 1.05
          repeat {
          ddlminimum <- departnoyau(dfmini,x,kernel,contr.sp$dftobwitmax,n,p,dfobjectif)
          if (ddlminimum<0) break else dfmini <- dfmini-0.02
        }
          if (abs(ddlminimum-dfobjectif)>0.1) {
            ddlmaximum <- 1
            dfmaxi <- 1.2
            repeat {
              ddlmaximum <- departnoyau(dfmaxi,x,kernel,contr.sp$dftobwitmax,n,p,dfobjectif)
              if (ddlmaximum>0) break else dfmaxi <- dfmaxi+0.2
            }
            if (abs(ddlmaximum-dfobjectif)>0.1) {
              res <- uniroot(departnoyau,c(dfmini,dfmaxi),tol=contr.sp$accuracy,x=x,kernel=kernel,dftobwitmax=contr.sp$dftobwitmax,n=n,p=p,dfobjectif)
              df <- res$root
            } else df <- dfmaxi
          } else df <- dfmini
        }
        bandwidth <- bwchoice(x,df,kernel,contr.sp$dftobwitmax)
        listeA <- calcA(X=x,bx=bandwidth,kernelx=kernel)
      }
    }
    dfstart <- listeA$df
    listeA.eig <- eigen(listeA$A,symmetric=TRUE)
    tPADmdemiY <- t(listeA.eig$vectors*(1/listeA$Ddemi))%*%y
    eigenvaluesA <- listeA.eig$values
    if (any(zapsmall(eigenvaluesA-1,digits=9)==0)) {
      ddlmini <-  sum(zapsmall(eigenvaluesA-1,digits=9)==0)
    } else ddlmini <- NA
    if (any(zapsmall(eigenvaluesA,digits=9)==0)) {
      index0 <-  which(zapsmall(eigenvaluesA,digits=9)==0)[1]
    } else index0 <- NA
    DdemiPA <- (listeA$Ddemi*listeA.eig$vectors)
    rm(listeA.eig)
    if (is.null(iter)) {
      if (Kmax<=Kmin) stop("Kmax should be greater than Kmin\n")
      if (any(c(Kmin,Kmax)<=0)) stop("Kmin and Kmax should be greater than 0")
      if (criterion%in%c("rmse","map")) {
        if (cv$bwchange) {
          if (is.null(df)) stop("df needs to be set\n")
          bx <- NULL } else bx <- bandwidth
        if (contr.sp$exhaustive) {
          choixk <- iterchoiceAcve(x,y,bx,df,kernel,ddlmini,cv$ntest,cv$ntrain,cv$Kfold,cv$type,cv$npermut,cv$seed,Kmin,Kmax)
          allcrit <- switch(criterion,rmse=choixk$rmse,map=choixk$map)
          iter <- (Kmin:Kmax)[which.min(allcrit)]
        } else {
          prov <- iterchoiceAcv(x,y,bx,df,kernel,ddlmini,cv$ntest,cv$ntrain,cv$Kfold,cv$type,cv$npermut,cv$seed,Kmin,Kmax,criterion,contr.sp$fraction)
          iter <- prov$iter
          choixk <- prov$objective
          if (criterion=="rmse") choixk <- sqrt(choixk)
        }
      } else {
        if (contr.sp$exhaustive) {
          choixk <- iterchoiceAe(y,Kmin:Kmax,eigenvaluesA,tPADmdemiY,DdemiPA,ddlmini,contr.sp$dfmaxi)
          allcrit <- switch(criterion,aic=choixk$aic,aicc=choixk$aicc,gcv=choixk$gcv,bic=choixk$bic,gmdl=choixk$gmdl)
          iter <- (Kmin:Kmax)[which.min(allcrit)]
        } else {
          prov <- iterchoiceA(n,Kmin,Kmax,eigenvaluesA,tPADmdemiY,DdemiPA,ddlmini,contr.sp$dfmaxi,y,criterion,contr.sp$fraction)
          iter <- prov$iter
          choixk <- prov$objective
        }
      }
      if ((((Kmax-iter)/(Kmax-Kmin)<1e-5)|(Kmax-iter)<3)&(!contr.sp$exhaustive)) warning(paste("Number of iterations is chosen close to the boundary of grid search: ",Kmax,".\nIncrease the maximum number of iterations or use contr.sp$exhaustive search\n",sep=""))
      if ((iter==max(Kmax))&(contr.sp$exhaustive)) warning(paste("Number of iterations is chosen at the boundary of grid search: ",Kmax,".\nIncrease the maximum number of iterations\n",sep=""))
    } else {
      criterion="user"
      choixk <- NULL
    }
    beta <- betaA(n,eigenvaluesA,tPADmdemiY,DdemiPA,ddlmini,k=iter,index0)
    listefit <- fittedA(n,eigenvaluesA,tPADmdemiY,DdemiPA,ddlmini,k=iter)
    residuals <- y-listefit$fit
  }
  if (smoother=="tps") {
    bandwidth <- contr.sp$bandwidth
    if (length(df)>1) stop("only one df is possible with Thin Plate Splines\n")
    ddlmini <- choose(contr.sp$m+p-1,contr.sp$m-1)
    if (is.null(bandwidth)) {
      lambda <- lambdachoice(x,ddlmini*df,m=contr.sp$m,itermax=contr.sp$dftobwitmax)
    } else {
      lambda <- bandwidth
    }
    S1 <- tpssmoother(x, y,lambda=lambda,m=contr.sp$m)
    vp1.S1 <- eigen(S1$H,symmetric=TRUE)
    U <- vp1.S1$vect
    eigenvaluesS1 <- vp1.S1$values
    rm(vp1.S1)
    tUy <- as.vector(crossprod(U, y))
    if (any(zapsmall(eigenvaluesS1,digits=9)==0)) {
      index0 <-  which(zapsmall(eigenvaluesS1,digits=9)==0)[1]
    } else index0 <- NA
    dfstart <- sum(eigenvaluesS1)
    if (is.null(iter)) {
      if (Kmax<=Kmin) stop("Kmax should be greater than Kmin\n")
      if (any(c(Kmin,Kmax)<=0)) stop("Kmin and Kmax should be greater than 0")
      if (criterion%in%c("rmse","map")) {
        if (cv$bwchange) {
          if (is.null(df)) stop("df needs to be set\n")
          lambda <- NULL }
        if (contr.sp$exhaustive) {
          choixk <- iterchoiceS1cve(x,y,lambda,df,ddlmini,cv$ntest,cv$ntrain,cv$Kfold,cv$type,cv$npermut,cv$seed,Kmin,Kmax,contr.sp$m)
          allcrit <- switch(criterion,rmse=choixk$rmse,map=choixk$map)
          iter <- (Kmin:Kmax)[which.min(allcrit)]
        } else {
          prov <- iterchoiceS1cv(x,y,lambda,df,ddlmini,cv$ntest,cv$ntrain,cv$Kfold,cv$type,cv$npermut,cv$seed,Kmin,Kmax,criterion,contr.sp$m,contr.sp$fraction)
          iter <- prov$iter
          choixk <- prov$objective
          if (criterion=="rmse") choixk <- sqrt(choixk)
        }
      } else {
        if (contr.sp$exhaustive) {
          choixk <- iterchoiceS1e(y,Kmin:Kmax,tUy,eigenvaluesS1,ddlmini,contr.sp$dfmaxi) 
          allcrit <- switch(criterion,aic=choixk$aic,aicc=choixk$aicc,gcv=choixk$gcv,bic=choixk$bic,gmdl=choixk$gmdl)
          iter <- (Kmin:Kmax)[which.min(allcrit)]
        } else {
          prov <- iterchoiceS1(n,Kmin,Kmax,tUy,eigenvaluesS1,ddlmini,contr.sp$dfmaxi,y,criterion,contr.sp$fraction)
          iter <- prov$iter
          choixk <- prov$objective
        }
      }
      if ((((Kmax-iter)/(Kmax-Kmin)<1e-5)|(Kmax-iter)<3)&(!contr.sp$exhaustive)) warning(paste("Number of iterations is chosen close to the boundary of grid search: ",Kmax,".\nIncrease the maximum number of iterations or use contr.sp$exhaustive search\n",sep=""))
      if ((iter==max(Kmax))&(contr.sp$exhaustive)) warning(paste("Number of iterations is chosen at the boundary of grid search: ",Kmax,".\nIncrease the maximum number of iterations\n",sep=""))
    } else {
      criterion="user"
      choixk <- NULL
    }
      listebeta <- betaS1(n,U,tUy,eigenvaluesS1,ddlmini,iter,lambda,S1$Sgu,S1$Qgu,index0)
    beta <- list(d=listebeta$dgub,c=listebeta$cgub)
    listefit <- fittedS1(n,U,tUy,eigenvaluesS1,ddlmini,iter)
    residuals <- y- listefit$fit
  }
  res <- list(beta=beta,residuals=residuals,fitted=listefit$fit,iter=iter,initialdf=dfstart,finaldf=listefit$trace,bandwidth=bandwidth,call=list(x=x,y=y,criterion=criterion,smoother=smoother,kernel=kernel,p=p,m=contr.sp$m),criteria=choixk)    
  class(res) <- c("ibr", "list")
  return(res)
}
