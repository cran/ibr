kernelKx <- function(X,Xetoile=NULL,bx,kernelx="g"){
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  namekernel <- c("g","e","q","u")
  kernelint <- which(namekernel==kernelx)
  if (length(kernelint)!=1) stop(paste("kernelx must be in",paste(namekernel, collapse=", ")))
  if (is.null(Xetoile)) {
    K <- .Call(ibr_Kmatrix, X, 0, bx, as.integer(c(n, p, 0, kernelint, 1)))
  } else {
    netoile <- nrow(Xetoile)
    K <- .Call(ibr_Kmatrix,  X, as.matrix(Xetoile), bx, as.integer(c(n, p, netoile, kernelint, 0)))
  }
 return(K)
}
