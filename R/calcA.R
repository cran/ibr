calcA <- function(X,bx,kernelx="g") {
  K <- kernelKx(X, NULL, bx, kernelx)
  listeS <- .Call(ibr_Amatrix, K)
  listeS[[3]] <- sum(diag(K)/listeS[[2]])
  listeS[[2]] <- 1/sqrt(listeS[[2]])
  names(listeS) <- c("S", "Ddemi", "df")
  return(listeS)
}

