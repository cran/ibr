#ifndef DSSPLINES_H_
#define DSSPLINES_H_
#include <R.h>
#include <Rinternals.h>
SEXP semikerlog(SEXP rx, SEXP rxk, SEXP rksurdeux, SEXP rparamsint);
SEXP semikerpow(SEXP rx, SEXP rxk, SEXP rksurdeux, SEXP rparamsint);
SEXP polynom(SEXP rparamsint, SEXP rdes);
#endif // DSSPLINES_H_
