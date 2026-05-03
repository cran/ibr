kernelSx <- function(X,Xetoile=NULL,bx,kernelx="g"){
  K <- kernelKx(X,Xetoile,bx,kernelx)
  S <- .Call(ibr_Smatrix, K)
  return(S)
}
