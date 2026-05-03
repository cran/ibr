#ifndef PRODUCT_H_
#define PRODUCT_H_
#include <R.h>
#include <Rinternals.h>
SEXP productandcrit(SEXP rS, SEXP ry, SEXP rcrits, SEXP rDm1, SEXP rparamsint); 
SEXP product(SEXP rS, SEXP ry, SEXP rparamsint);
#endif // PRODUCT_H_
