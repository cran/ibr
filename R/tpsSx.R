tpsSx <- function(X,Xetoile,lambda,m=2) {
  n <- nrow(X)
  d <- ncol(X)
  p <- 2*m-d
  n <- nrow(X)
  netoile <- nrow(Xetoile)
  d <- ncol(X)
  if (ncol(Xetoile)!=d) stop("the number of variables do not match")
  Sgu <-  fields.mkpoly(Xetoile, m = m)
  Kgu <- Rad.cov(Xetoile,X,p=p)
  return(list(Sgu=Sgu,Qgu=Kgu))
}
