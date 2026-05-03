departnoyau <- function(df,x,kernel,dftobwitmax,n,p,dfobjectif) {
  bx <- bwchoice(x,df,kernel,dftobwitmax)
  H <- kernelSx(x, NULL, bx, kernel)
  trace <- sum(diag(H))
  return(trace)
}
