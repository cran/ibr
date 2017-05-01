#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void gaustotal(double *ax, double *bx, double *x, int *nx, int *px, double *Tol, int *Maxit, double *objectif, int *dftotal, double *K, double *Ddemi, double *dfstart, double *bandwidth );
extern void polynom(int *m, int *n, int *dim, double *des, int *lddes, int *npoly, double *t, int *ldt, int *wptr, int *info);
extern void semikerlog(double *x, double *xk, int *nx, int *nxk,double *ksurdeux, int *px, int *negatif, double *dista, int *symmetric);
extern void semikerpow(double *x, double *xk, int *nx, int *nxk,double *ksurdeux, int *px, int *negatif, double *dista, int *symmetric);

extern void regg(double *x, int *nx,double *y, double *bw, double *valx, int *nvalx, double *regx, double *df);
extern void rege(double *x, int *nx,double *y, double *bw, double *valx, int *nvalx, double *regx, double *df);
extern void regq(double *x, int *nx,double *y, double *bw, double *valx, int *nvalx, double *regx, double *df);

extern void reggcv(double *x, int *nx,double *y, double *bw, int *nbw, int *effold, int *neffold, double *sse, double *sap );
extern void regecv(double *x, int *nx,double *y, double *bw, int *nbw, int *effold, int *neffold, double *sse, double *sap );
extern void regqcv(double *x, int *nx,double *y, double *bw, int *nbw, int *effold, int *neffold, double *sse, double *sap );

extern void regpolg(double *x, int *nx,double *y, double *bw, double *valx, int *nvalx, double *regx,  double *df, double *deriv);
extern void regpolq(double *x, int *nx,double *y, double *bw, double *valx, int *nvalx, double *regx,  double *df, double *deriv);
extern void regpole(double *x, int *nx,double *y, double *bw, double *valx, int *nvalx, double *regx,  double *df, double *deriv);
  
extern void regpolgcv(double *x, int *nx,double *y, double *bw, int *nbw, double *effold, int *neffold, double *sse, double *sap );
extern void regpolecv(double *x, int *nx,double *y, double *bw, int *nbw, double *effold, int *neffold, double *sse, double *sap );
extern void regpolqcv(double *x, int *nx,double *y, double *bw, int *nbw, double *effold, int *neffold, double *sse, double *sap );

  

static const R_CMethodDef CEntries[] = {
    {"gaustotal",  (DL_FUNC) &gaustotal,  13},
    {"polynom",    (DL_FUNC) &polynom,    12},
    {"semikerlog", (DL_FUNC) &semikerlog,  9},
    {"semikerpow", (DL_FUNC) &semikerpow,  9},
    {"regg", (DL_FUNC) &regg,  8},
    {"regq", (DL_FUNC) &regq,  8},
    {"rege", (DL_FUNC) &rege,  8},
    {"reggcv", (DL_FUNC) &reggcv,  9},
    {"regqcv", (DL_FUNC) &regqcv,  9},
    {"regecv", (DL_FUNC) &regecv,  9},
    {"regpolg", (DL_FUNC) &regpolg,  9},
    {"regpolq", (DL_FUNC) &regpolq,  9},
    {"regpole", (DL_FUNC) &regpole,  9},
    {"regpolgcv", (DL_FUNC) &regpolgcv,  9},
    {"regpolqcv", (DL_FUNC) &regpolqcv,  9},
    {"regpolecv", (DL_FUNC) &regpolecv,  9},
    {NULL, NULL, 0}
};

void R_init_ibr(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
