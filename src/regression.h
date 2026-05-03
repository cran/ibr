#ifndef REGRESSION_H_
#define REGRESSION_H_
#include <R.h>
#include <Rinternals.h>
SEXP npreg(SEXP rx, SEXP ry, SEXP rvalx, SEXP rbw, SEXP rparamsint);
SEXP npregcv(SEXP rx, SEXP ry, SEXP rbw, SEXP reffold, SEXP rparamsint);
#endif // REGRESSION_H_
