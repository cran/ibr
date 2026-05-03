ibr.fit <- function(x,y,criterion="gcv",df=1.5,Kmin=1,Kmax=1e+06,smoother="k",kernel="g",method="eigen",rank=NULL,control.par=list(),cv.options=list()) {
  cl <- match.call() 
  crit <- c("aic","aicc","gcv","bic","gmdl","rmse","map")
  crite <- pmatch(criterion,crit)
  if (any(is.na(crite))) stop(paste("parameter criterion must be in",paste(crit,collapse=", "))) else criterion <- crit[crite]
  methode <- c("eigen", "product")
  nummethod <- pmatch(method, methode)
  if (any(is.na(nummethod))) stop(paste("parameter method must be in",paste(methode,collapse=", "))) else method <- methode[nummethod]
  if ((any(crite>5))&&(length(crite)>2)) stop("RMSE or MAP must be used without any other criterion")
  iterautre <- NULL
  choixkautre <- NULL
  ##  criterion <- match.arg(criterion,crit)
  smoothertable <- c("k","tps","ds","lrtps","lrds")
  smoother <- match.arg(smoother,smoothertable)
  lowrank <- ifelse(substr(smoother,0,2)=="lr",TRUE,FALSE)
  if ((lowrank)&(method!="eigen"))  {
    warning("For lowrank splines only \"eigen\" method (one need to diagonalize to get smoother)")
    method <- "eigen"
  }
  if (!is.matrix(x)) {
    x <- data.matrix(x)
    warning("x is coerced to matrix by data.matrix function\n")
  }
  n <- nrow(x)
  p <- ncol(x)
  if (length(y)!=n) stop("number of observations in x and y do not match\n")
  if (smoother=="k") contr.sp <- list(bandwidth=NULL,iter=NULL,really.big=FALSE,
                                      dftobwitmax=1000,exhaustive=FALSE,m=NULL,s=NULL,dftotal=FALSE,
                                      accuracy=0.01,dfmaxi=2*n/3,fraction=c(100, 200, 500, 1000, 5000,10^4,5e+04,1e+05,5e+05,1e+06),
                                      scale=FALSE,aggregcrit="no",aggregfun=function(x) {floor(stats::median(x[criterion]))})
  else  contr.sp <- list(bandwidth=NULL,iter=NULL,really.big=FALSE,
                         dftobwitmax=1000,exhaustive=FALSE,m=NULL,s=NULL,dftotal=FALSE,
                         accuracy=0.01,dfmaxi=2*n/3,fraction=c(100, 200, 500, 1000, 5000,10^4,5e+04,1e+05,5e+05,1e+06),
                         scale=TRUE,aggregcrit="no",aggregfun=function(x) {floor(stats::median(x[criterion]))})
  contr.sp[(names(control.par))] <- control.par
  if (!(is.logical(contr.sp$dftotal))) stop("contr.sp$dftotal must be logical\n")
  if ((!is.null(contr.sp$bandwidth))&(!is.numeric(contr.sp$bandwidth) || any(contr.sp$bandwidth<0))) stop("invalid bandwidth\n")
  if ((!is.null(contr.sp$iter))&(!is.numeric(contr.sp$iter) || (contr.sp$iter<0) || (floor(contr.sp$iter)!=contr.sp$iter))) stop("invalid number of iterations\n")
  iter <- contr.sp$iter
  if ((contr.sp$dfmaxi<=0)|(contr.sp$dfmaxi>n)) stop("invalid dfmaxi\n")
  crit <-c("recalc","no","aggregation")
  crite <- pmatch(contr.sp$aggregcrit,crit)
  if (is.na(crite)) stop(paste("control.par$criterion must be in",crit))
  if (length(criterion)==2&&(all(criterion%in%c("map","rmse")))&&(!contr.sp$exhaustive)) stop("when RMSE and MAP are both selected, exhaustive search (in control.par list) must be chosen")
  if ((smoother=="tps")|(smoother=="lrtps")) {
    contr.sp$s <- 0
    if (is.null(contr.sp$m)) contr.sp$m <- floor(p/2)+1 else {
                                                          if (contr.sp$m<=(p/2)) stop("order of thin plate splines is invalid (need to be greater than p/2)\n")
                                                        }
    if (!is.numeric(contr.sp$m) || (contr.sp$m<0) || (floor(contr.sp$m)!=contr.sp$m)) stop("invalid spline order\n")
    ddlmin <- choose(contr.sp$m+p-1,contr.sp$m-1)
    if (ddlmin>=n) stop(paste("ddl min is equal to",ddlmin,"and the number observations is",n))
  }
  if (smoother=="k") {
    contr.sp$m <- NULL
  }
  smoothobject <- NULL
  if ((smoother=="ds")|(smoother=="lrds")) {
    if (is.null(contr.sp$m)) contr.sp$m <- 2 ## default penalty order 2
    if (is.null(contr.sp$s)) contr.sp$s <- (p-1)/2 ## default pseudo cubic
    if ((!is.numeric(contr.sp$m))|(!is.numeric(contr.sp$s))) stop("contr.par$m or contr.par$s is not numeric...\n")
    contr.sp$m <- round(contr.sp$m)     ## m is integer
    contr.sp$s <- round(contr.sp$s*2)/2 ## s is in halfs
    if (contr.sp$m< 1) contr.sp$m <- 1  ## m > 0
    ## check that -p/2 < s < p/2...
    if (contr.sp$s >= p/2) { 
      contr.sp$s <- (p-1)/2
      warning("contr.par$s value reduced")
    } 
    if (contr.sp$s <= -p/2) { 
      contr.sp$s <- -(p-1)/2
      warning("contr.par$s value increased")
    }
    
    ## m + s > p/2 for continuity...
    if ((contr.sp$m+contr.sp$s)<=p/2) {
      contr.sp$s <- 1/2 + p/2 - contr.sp$m
      if (contr.sp$s>=p/2) stop("No suitable contr.par$s try increasing contr.par$m")
      warning("contr.par$s value modified to give continuous function")
    }
    ddlmin <- choose(contr.sp$m+p-1,contr.sp$m-1)
    if (ddlmin>=n) stop(paste("ddl min is equal to",ddlmin,"and the number observations is",n))
  }
  
  if (lowrank) {
    if (is.null(rank)) stop("rank argument for lowrank splines must be chosen...")
    if (rank>n) stop(paste("rank argument must be less than",n))
    bs <- substr(smoother,3,4)
    listvarx <- colnames(x)
  } 
  moy <- NULL
  ec <- NULL
  if (contr.sp$scale) {
    if (smoother=="k") warning("when using kernel smoother, you do not need to scale\n")
    if (lowrank) ec <- apply(x,2, sd)*sqrt((n-1)/n) else  ec <- apply(x,2, sd)
    x <- scale(x,scale=ec)
    moy <- attr(x,"scaled:center")
  }
  if (all(criterion%in%c("rmse","map"))) {
    cv <- list(bwchange=FALSE,ntest=floor(nrow(x)/10),ntrain=NULL,Kfold=FALSE,type="random",seed=NULL,npermut=20)
    cv[(names(cv.options))] <- cv.options
    if (!all(sapply(cv[1],is.logical))) stop("invalid cv$bwchange or cv$Kfold: must be logical\n")
    if (!all(sapply(cv[c(2,3,6,8)], FUN=function(x) is.numeric(x)||is.null(x)))) stop("invalid cv parameters: must be numeric or NULL\n")
    if (any(names(cv.options)=="ntrain")) cv$ntest <- NULL
  } else cv <- NULL
  if (!(lowrank)&&((n>1000)&(! contr.sp$really.big))) stop("number of observations is greater than 1000, set control.par$really.big to TRUE if you really want to do the requested calculations (but computational time -eigen decomposition- could be prohibitive)\n")
  if (smoother=="k") {
    kern <- c("g", "e", "q", "u")
    kernel <- match.arg(kernel,kern)
    kernelint <- which(kernel==kern)
    m <- NULL
    if (!is.null(contr.sp$bandwidth)) {
      if (length(contr.sp$bandwidth)==1) {
        bandwidth <- rep(contr.sp$bandwidth,p)
      } else {
        if (length(contr.sp$bandwidth)!=p)  stop(paste("the length of bandwidth vector have to be",p,"or 1\n"))
        bandwidth <- contr.sp$bandwidth
      }
    } else {
      if (df<=1) stop("degree of freedom should be greater than 1\n")   
      departbw <- apply(x,2,FUN=function(z) 3*abs(diff(range(z))))
      bandwidth <- .Call(ibr_choosebw,
                         as.double(departbw),
                         if (contr.sp$dftotal) as.double(rep(1e-10,p)) else
                                                                         as.double(1e-10),
                         as.double(x),
                         c(as.double(.Machine$double.eps^0.25),
                           as.double(df)),
                         as.integer(c(n, p, contr.sp$dftobwitmax,
                                      contr.sp$dftotal, kernelint))
                         )
    }
    K <- .Call(ibr_Kmatrix, x, 0, bandwidth, as.integer(c(n, p, n, kernelint, 1)))
    listeS <- .Call(ibr_Amatrix, K)
    listeS[[3]] <- sum(diag(K)/listeS[[2]])
    listeS[[2]] <- 1/sqrt(listeS[[2]])
    names(listeS) <- c("S", "Ddemi", "df")
    rm(K)
    dfstart <- listeS$df
  } # end kernel
  if ((smoother=="tps")|(smoother=="ds")|lowrank) {
    bandwidth <- contr.sp$bandwidth
    if (length(df)>1) stop("only one df is possible with Splines\n")
    ddlmini <- choose(contr.sp$m+p-1,contr.sp$m-1)
    if (lowrank) {
      if (is.null(bandwidth)) {
        lambda <- lambdachoicelr(x,ddlmini*df,m=contr.sp$m,contr.sp$s,rank,itermax=contr.sp$dftobwitmax,bs,listvarx)
        bandwidth <- lambda
      } else {
        lambda <- bandwidth
      }
      S2 <- lrsmoother(x,bs,listvarx,lambda=lambda,m=contr.sp$m,s=contr.sp$s,rank)
    } else {
      if (is.null(bandwidth)) {
        lambda <- lambdachoice(x,ddlmini*df,m=contr.sp$m,contr.sp$s,itermax=contr.sp$dftobwitmax,smoother)
        bandwidth <- lambda
      } else {
        lambda <- bandwidth
      }
      listeS <- dssmoother(x, y,lambda=lambda,m=contr.sp$m,s=contr.sp$s)
      dfstart <- sum(diag(listeS$S))
    }
  }
  if (method=="eigen") {
    if (lowrank) {
      ## diagonalization already done
      eigenvaluesS <- S2$values
      dfstart <- sum(eigenvaluesS)
      ## almost 0 eigenvalues
      if (any(zapsmall(eigenvaluesS,digits=9)==0)) {
        index0 <-  which(zapsmall(eigenvaluesS,digits=9)==0)[1]
      } else index0 <- rank+1
      eigenvaluesS <- c(eigenvaluesS,rep(0,n-rank))
      eigenvaluesS[eigenvaluesS<0] <- 0
      U <- S2$vectors
      Rm1U <- S2$Rm1U
      tUy <- c(as.vector(crossprod(U, y)),rep(0,n-rank))
      smoothobject <- S2$smoothobject
      rm(S2)
    } else {
      ## diagonalization
      listeS.eig <- eigen(listeS$S,symmetric=TRUE)
      eigenvaluesS <- listeS.eig$values
      if (any(eigenvaluesS<(-1e-10))) stop("Some eigenvalues of the Kernel smoother matrix are negative, it will explode")
      if (smoother=="k") {
        tPADmdemiY <- t(listeS.eig$vectors*(1/listeS$Ddemi))%*%y
        DdemiPA <- (listeS$Ddemi*listeS.eig$vectors)
      } else {
        U <- listeS.eig$vectors
        tUy <- as.vector(crossprod(U, y))
      }
      ## ddlmini for kernel
      if (smoother=="k"){
        if (any(zapsmall(eigenvaluesS-1,digits=9)==0)) {
          ddlmini <-  sum(zapsmall(eigenvaluesS-1,digits=9)==0)
        } else ddlmini <- max(eigenvaluesS)
      }
      ## almost 0 eigenvalues
      if (any(zapsmall(eigenvaluesS,digits=9)==0)) {
        index0 <-  which(zapsmall(eigenvaluesS,digits=9)==0)[1]
      } else index0 <- NA   
      rm(listeS.eig)  
    }
    if (is.null(iter)) {
      if (Kmax<=Kmin) stop("Kmax hould be greater than Kmin\n")
      if (any(c(Kmin,Kmax)<=0)) stop("Kmin and Kmax should be greater than 0")
      if (all(criterion%in%c("rmse","map"))) {
        ## iter k chosen by CV
        if (cv$bwchange) {
          if (is.null(df)) stop("df needs to be set\n")
          if (smoother=="k") bx <- NULL else lambda <- NULL
        } else {
          if (smoother=="k") bx <- bandwidth
        }
        if (contr.sp$exhaustive) {
          ## exhaustive
          if (smoother=="k") {
            choixkautre <- iterchoiceAcve(x,y,bx,df,kernel,ddlmini,cv$ntest,
                                          cv$ntrain,cv$Kfold,cv$type,cv$npermut,
                                          cv$seed,Kmin,Kmax)
          } else {
            if (lowrank)     {
              choixkautre <- iterchoiceS1lrcve(x,y,lambda,rank,bs,listvarx,df,ddlmini,cv$ntest,cv$ntrain,cv$Kfold,cv$type,cv$npermut,cv$seed,Kmin,Kmax,contr.sp$m,contr.sp$s)
            } else {
              choixkautre <- iterchoiceS1cve(x,y,lambda,df,ddlmini,cv$ntest,cv$ntrain,cv$Kfold,cv$type,cv$npermut,cv$seed,Kmin,Kmax,contr.sp$m,contr.sp$s)
            }
          }
          iterautre <-  (Kmin:Kmax)[unlist(lapply(choixkautre[1:5],FUN=which.min))]
          names(iterautre) <-  names(choixkautre[1:5])
          iter <- switch(contr.sp$aggregcrit,
                         no=iterautre[criterion[1]],
                         aggregation=contr.sp$aggregfun(iterautre),
                         recalc= { Kmax2 <- stats::median(10^(floor(log10(iterautre))+2)) ;
                           iterautre[1] <- (Kmin:Kmax2)[which.min(choixkautre[[criterion[1]]])] ;
                           iterautre[1] })
          if (contr.sp$aggregcrit=="aggregation") {
            choixk <- NA
            names(choixk) <- "aggregation"
          } else {
            choixk <- switch(criterion[1],aic=choixkautre$aic[iter],aicc=choixkautre$aicc[iter],gcv=choixkautre$gcv[iter],bic=choixkautre$bic[iter],gmdl=choixkautre$gmdl[iter])
            names(choixk) <- criterion[1]
          }
        } else {
          ## optimize 
          if (smoother=="k") {
            prov <- iterchoiceAcv(x,y,bx,df,kernel,ddlmini,cv$ntest,cv$ntrain,cv$Kfold,cv$type,cv$npermut,cv$seed,Kmin,Kmax,criterion,contr.sp$fraction)
          } else {
            if (lowrank)     {
              prov <- iterchoiceS1lrcv(x,y,lambda,rank,bs,listvarx,df,ddlmini,cv$ntest,cv$ntrain,cv$Kfold,cv$type,cv$npermut,cv$seed,Kmin,Kmax,criterion,contr.sp$m,contr.sp$s,contr.sp$fraction)
            } else {
              prov <- iterchoiceS1cv(x,y,lambda,df,ddlmini,cv$ntest,cv$ntrain,cv$Kfold,cv$type,cv$npermut,cv$seed,Kmin,Kmax,criterion,contr.sp$m,contr.sp$s,contr.sp$fraction)
            }
          }
          iter <- prov$iter
          choixkautre <- prov$objective
          names(choixk) <- criterion[1]
          if (criterion=="rmse") choixk <- sqrt(choixk)
        }
      } else {
        ## iter k chosen by AIC/BIC etc
        if (contr.sp$exhaustive) {
          ## exhaustive
          if (smoother=="k") {
            choixkautre <- iterchoiceAe(y,Kmin:Kmax,eigenvaluesS,tPADmdemiY,DdemiPA,ddlmini,contr.sp$dfmaxi)
            
          } else {
            choixkautre <- iterchoiceS1e(y,Kmin:Kmax,tUy,eigenvaluesS,ddlmini,contr.sp$dfmaxi)
          }
          iterautre <- (Kmin:Kmax)[unlist(lapply(choixkautre[c("aic","aicc","gcv","bic","gmdl")],FUN=which.min))]
          names(iterautre) <-  names(choixkautre[1:5])
          iter <- switch(contr.sp$aggregcrit,
                         no=iterautre[ criterion[1] ],
                         aggregation=contr.sp$aggregfun(iterautre),
                         recalc= { Kmax2 <- stats::median(10^(floor(log10(iterautre))+2)) ;
                           iterautre[criterion[1]] <- (Kmin:Kmax2)[ which.min(choixkautre[[ criterion[1] ]])] ;
                           iterautre[criterion[1]] })
          if (contr.sp$aggregcrit=="aggregation") {
            choixk <- NA
            names(choixk) <- "aggregation"
          } else {
            choixk <- switch(criterion[1],aic=choixkautre$aic[iter],aicc=choixkautre$aicc[iter],gcv=choixkautre$gcv[iter],bic=choixkautre$bic[iter],gmdl=choixkautre$gmdl[iter])
          }
          ## 
        } else {
          ## optimize
          if (smoother=="k") {
            prov <- iterchoiceA(n,Kmin,Kmax,eigenvaluesS,tPADmdemiY,DdemiPA,ddlmini,contr.sp$dfmaxi,y,criterion[1],contr.sp$fraction)
          }  else {
            prov <- iterchoiceS1(n,Kmin,Kmax,tUy,eigenvaluesS,ddlmini,contr.sp$dfmaxi,y,criterion[1],contr.sp$fraction)
          }
          iter <- prov$iter
          choixk <- prov$objective
          names(choixk) <- criterion[1]
          iterautre <- iter
          choixkautre <- prov$objective
          if (length(criterion)>1) {
            ## several criteria
            for (i in 2:length(criterion)) {
              if (smoother=="k") {
                prov <- iterchoiceA(n,Kmin,Kmax,eigenvaluesS,tPADmdemiY,DdemiPA,ddlmini,contr.sp$dfmaxi,y,criterion[i],contr.sp$fraction)
              } else {
                prov <- iterchoiceS1(n,Kmin,Kmax,tUy,eigenvaluesS,ddlmini,contr.sp$dfmaxi,y,criterion[i],contr.sp$fraction)
              }
              iterautre <- c(iterautre,prov$iter)
              choixkautre <- c(choixkautre,prov$objective)
            }
            names(iterautre) <- criterion
            iter <- switch(contr.sp$aggregcrit,no=iterautre[1],
                           aggregation={
                             choixk <-  "aggregation"
                             contr.sp$aggregfun(iterautre) }, recalc= {
                               Kmax2 <- stats::median(10^(floor(log10(iterautre))+2))
                               if (smoother=="k") {
                                 prov <- iterchoiceA(n,Kmin,Kmax2,eigenvaluesS,tPADmdemiY,DdemiPA,ddlmini,contr.sp$dfmaxi,y,criterion[1],contr.sp$fraction)
                               } else {
                                 prov <- iterchoiceS1(n,Kmin,Kmax2,tUy,eigenvaluesS,ddlmini,contr.sp$dfmaxi,y,criterion[1],contr.sp$fraction)
                               }
                               choixk <- prov$objective
                               iterautre[1] <- prov$iter
                               prov$iter }) 
          } ## several criteria
        } ## exhaustive or optimize
      } ## iter K chosen by CV or Not
      if ((((Kmax-iter)/(Kmax-Kmin)<1e-5)|(Kmax-iter)<3)&(!contr.sp$exhaustive)) warning(paste("Number of iterations is chosen close to the boundary of grid search: ",Kmax,".\n  Increase the maximum number of iterations or use contr.sp$exhaustive search\n",sep=""))
      if ((iter==max(Kmax))&(contr.sp$exhaustive)) warning(paste("Number of iterations is chosen at the boundary of grid search: ",Kmax,".\n  Increase the maximum number of iterations\n",sep=""))
    } else {
      criterion <- "user"
      choixk <- NULL
    } ## iter K chosen by crit or user
  } else {
    ## product
    if (Kmin!=1) {
      Kmin <- 1
      warning("Kmin must be one for product: Kmin is set to 1 and proceed")
    }
    if (smoother=="k") {
      Ddemiy <- y/listeS[[2]]
    } else {
      Ddemiy <- y
    }
    if (!is.null(iter)) {
      ## prod until iter
      res <- .Call(ibr_product, listeS$S, as.double(Ddemiy), as.integer(iter))
      criterion <- "user"
      choixk <- NULL
      choixkautre <- NULL
      beta <- res[[1]]
      finaldf <- res[[2]]
    } else {
      if(all(criterion%in%c("aic","aicc","bic","gcv"))) {
        prov <- c("aic","aicc","bic","gcv") %in% criterion
        if (any(prov)) { 
          critsint <- which(prov)
          ## if (length(which(prov))>1) warning("only one criterion for method=\"product\"")
        } else {
          stop("Criterion must be in \"aic\",\"aicc\",\"gcv\",\"bic\" for \"product\" method")
        }
      } else {
        stop("criterion must be in aic, aicc, bic, gcv")
      }
      if (smoother=="k") {
        res <- .Call(ibr_productandcrit, listeS$S, as.double(Ddemiy),
                     as.integer(critsint),
                     as.double(listeS$Ddemi),
                     as.integer(c(Kmax, ifelse(smoother=="k", 0, 1))))
      } else {
        res <- .Call(ibr_productandcrit, listeS$S, as.double(Ddemiy),
                     as.integer(critsint),
                     as.double(0.0),
                     as.integer(c(Kmax, ifelse(smoother=="k", 0, 1))))
      }
      if (length(criterion)>1) { 
        myorder <- order(criterion)
        betaautre <- matrix(res[[1]], n, length(criterion))[, myorder, drop=FALSE]
        beta <-  betaautre[,1]
        iterautre <- res[[4]][myorder]
        names(iterautre) <- criterion
        choixkautre <- res[[3]][myorder]
        finaldf <-  res[[2]][myorder]
        if (contr.sp$aggregcrit=="no") {
          iter <- iterautre[1]
          myname <- criterion[1]
          choixk <- c(myname = choixkautre[1])
          finaldf <- finaldf[1]
        }
        if (contr.sp$aggregcrit=="aggregation") {
          iter <- contr.sp$aggregfun(iterautre)
          res <- .Call(ibr_product, listeS$S, as.double(Ddemiy), as.integer(iter))
          choixk <- c(aggregation = NA)
          finaldf <- res[[2]]
        }
        if (contr.sp$aggregcrit=="recalc") {
          Kmax2 <- stats::median(10^(floor(log10(iterautre))+2))
          prov <- (c("aic","aicc","bic","gcv")==criterion[1])
          critsint <- which(prov)
          res <- .Call(ibr_productandcrit, listeS$S, as.double(Ddemiy),
                       as.integer(critsint),
                       as.double(0.0),
                       as.integer(c(Kmax2, ifelse(smoother=="k", 0, 1))))
          beta <-  res[[1]]
          iterautre <- res[[4]]
          choixkautre <- res[[3]]
          names(choixkautre) <- criterion
          choixk <- choixkautre
          iter <- iterautre
          finaldf <- res[[2]]
        }
      } else {
        # one crit only
        beta <-  res[[1]]
        iterautre <- res[[4]]
        choixkautre <- res[[3]]
        names(choixkautre) <- criterion
        choixk <- choixkautre
        iter <- iterautre
        finaldf <- res[[2]]
      }
    }
  } # end of product
  ## beta
  if (method=="eigen") {
    if (smoother=="k") {
      beta <- betaA(n,eigenvaluesS,tPADmdemiY,DdemiPA,ddlmini,k=iter,index0)
      listefit <- fittedA(n,eigenvaluesS,tPADmdemiY,DdemiPA,ddlmini,k=iter)
    } else {
      if (lowrank) {
        beta <- betaS1lr(n,U,tUy,eigenvaluesS,ddlmini,iter,lambda,rank,Rm1U,index0)
        listefit <- fittedS1lr(n,U,tUy,eigenvaluesS,ddlmini,iter,rank)
        if (listefit$trace>=0.99*rank) warning(paste("Final rank is chosen equal to the lowrank of spline: ",rank,".\n  Increase the rank argument\n",sep=""))
      } else {
        listebeta <- betaS1(n,U,tUy,eigenvaluesS,ddlmini,iter,lambda,listeS$Sgu,listeS$Qgu,index0)
        beta <- list(d=listebeta$dgub,c=listebeta$cgub)
        listefit <- fittedS1(n,U,tUy,eigenvaluesS,ddlmini,iter)
      }
    }
    finaldf <- listefit$trace
    fitted <- listefit$fit
  } else {
    if (smoother=="k") {
      fitted <- matrix(rep(listeS$Ddemi, n)*listeS$S,n,n)%*%beta
      beta <- listeS$Ddemi*beta
      } else {
       fitted <- listeS$S%*%beta
      }
  }
  residuals <- y - fitted
  res <- list(beta=beta,residuals=residuals,fitted=fitted,iter=iter,initialdf=dfstart,
              finaldf=finaldf,bandwidth=bandwidth,call=cl,
              parcall=list(p=p,m=contr.sp$m,s=contr.sp$s,scaled=contr.sp$scale,mean=moy,sd=ec,critmethod=contr.sp$aggregcrit,method=method,
                           rank=rank,criterion=criterion,smoother = smoother, kernel = kernel,smoothobject=smoothobject,exhaustive=contr.sp$exhaustive),
              criteria=choixk,alliter=iterautre,allcriteria=choixkautre)    
  return(res)
}
