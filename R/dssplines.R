DuchonQ <- function(x,xk,m=2,s=0,symmetric=TRUE) {
  p <- ncol(x)
  nx <- nrow(x)
  if (!symmetric) {
    if (ncol(xk)!=p) stop("number of variables of x and xk must be the same\n")
    nxk <- nrow(xk)
  } else {
    nxk <- nx
    xk <- 0
  }
  k <- 2*m + 2*s - p
  negatif <- if ((1-2*((floor(k/2)+1)%%2))==-1) 1 else 0
  if (k%%2==0) {
    res <- .Call(ibr_semikerlog,as.double(x),as.double(xk),as.double(k/2),
                 as.integer(c(nx,p,nxk,negatif,symmetric)))
   } else {
    res <- .Call(ibr_semikerpow,as.double(x),as.double(xk),as.double(k/2),
                 as.integer(c(nx,p,nxk,negatif,symmetric)))
  }
 return(matrix(res,nrow=nx,ncol=nxk))
}

# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
DuchonS <- function(x, m = 2) {
    if (m < 1) 
        stop("'m' has to be larger than zero.")
    if (!is.matrix(x)) 
        x <- as.matrix(x)
    d <- ncol(x)
    n <- nrow(x)
    nterms<- choose((m + d -1),d)
    temp <- .Call(ibr_polynom, as.integer(c(m, n, d, n, nterms, n)),
        des = as.double(x))
    return(matrix(temp, nrow = n))
}
