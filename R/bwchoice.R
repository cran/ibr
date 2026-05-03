bwchoice <- function(X,objectif,kernelx="g",itermax=1000) {
  p <- ncol(as.matrix(X))
  res <- rep(0,p)
  if (length(objectif)==1) objectif <- rep(objectif,p)
  if (any(objectif<=1)) stop("degree of freedom should be greater than 1\n")
  choixddlparvar <- function(fenetre,X,objectif,kernel) {
    H <- kernelSx(X, NULL, fenetre, kernel)
    trace <- sum(diag(H))
    res <- trace-objectif
    return(res)
  }
  for (j in 1:ncol(X)) {
    depart <- 3*abs(diff(range(X[,j])))
    if (choixddlparvar(depart,X[,j],objectif[j],kernelx)>0) {
      while (choixddlparvar(depart,X[,j],objectif[j],kernelx)>0) {
        depart <- depart*2
      }
    }
    res[j] <- stats::uniroot(choixddlparvar,interval=c(depart,1e-10),X=X[,j],objectif=objectif[j],maxiter=itermax,kernel=kernelx)$root
  }
  return(res)
}
