#ifndef POLYNOMLOCAUX_H_
#define POLYNOMLOCAUX_H_
#include <R.h>
#include <Rinternals.h>
SEXP npregpol(SEXP rx, SEXP ry, SEXP rvalx, SEXP rbw, SEXP rparamsint);
SEXP npregpolcv(SEXP rx, SEXP ry, SEXP rbw, SEXP reffold, SEXP rparamsint);
#endif // POLYNOMLOCAUX_H_
