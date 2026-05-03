#ifndef SMOOTHER_H_
#define SMOOTHER_H_
#include <R.h>
#include <Rinternals.h>
SEXP Kmatrix(SEXP rx, SEXP rvalx, SEXP rbandwidth, SEXP rparamsint);
SEXP Smatrix(SEXP rK);
SEXP Amatrix(SEXP rK);
#endif // SMOOTHER_H_
