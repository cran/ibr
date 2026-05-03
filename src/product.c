#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h> // interrupt
#include <Rmath.h>
#include "product.h"

static R_INLINE double aic(double scr, double df, int n)
{
  double res;
  res = log(scr/n) +  2* df /n ;
  return res;
}
static R_INLINE double aicc(double scr, double df, int n)
{
  double res;
  res = log(scr/n) +  1 + 2*(df + 1)/(n - df -2) ;
  return res;
}

static R_INLINE double bic(double scr, double df, int n)
{
  double res;
  res = log(scr/n) +  log(n)* df /n ;
  return res;
}

static R_INLINE double gcv(double scr, double df, int n)
{
  double res;
  res = log(scr/n) - 2 * log(1 - df/n) ;
  return res;
}




static double calculcrit(double scr, double df, int n, int critnumber) {
  /* return the criterion */
  double res;
   switch(critnumber) {
  case 1:
    /* AIC */
    res= aic(scr, df, n);
    break;
  case 2:
    /* AICc */
    res= aicc(scr, df, n);
    break;
  case 3:
    /* BIC */
    res= bic(scr, df, n);
    break;
  case 4:
    /* GCV */
    res= gcv(scr, df, n);
    break;
   default:
     res=R_PosInf;
  }
  return res;
}


SEXP productandcrit(SEXP rS, /* Smoother matrix symmetric */
	                SEXP ry, /* data Y */
                    SEXP rcritsint, /* int vector of criteria to calculate */
                    SEXP rDmdemi, /* diag(D^{-1}) used if kernel */
                    SEXP rparamsint) 
{
  int k, i, j, n, *critsint, ncrits, niter, symmetric, *critsok, allok, lenImSk;
  double *S, *Dmdemi, *y, *gamma, *gammatemp, *gammacur, *ImSkold, *ImSk, *ImS, *critsold, *itercrits, *df;
  double *crits, dftemp, dfS=0.0, dfSold, scr, one=1.0, zero=0.0;
  const double MAXI=1e10;
  char *uplo="u", *side="r";
  SEXP rans, rgamma, rcrits, ritercrits, rdf;
  /* rparamsint
     - niter maximum number of iteration greater than 1
     - symmetric (1 if S symmetric, 0 if not)
   */
  niter = INTEGER(rparamsint)[0];
  symmetric = INTEGER(rparamsint)[1];
  ncrits=Rf_length(rcritsint);
  n=Rf_length(ry);
  /* ---------- output ------------  */
  /* rans: list  */
  /*   - rgamma: vector of double (size n * ncrits)
                 contains D^{1/2}\beta^{(k)}
                 k is chosen by each requested crit
       - rdf: vector of double (size ncrits), df at min of each crit
       - rcrits: vector of double (size ncrits), min of crit
       - ritercrits: vector of int (size ncrits), iter of min
   */
  rans=PROTECT(Rf_allocVector(VECSXP, 4));
  rgamma =  PROTECT(Rf_allocVector(REALSXP, n*ncrits));
  gamma =  REAL(rgamma);
  rcrits =  PROTECT(Rf_allocVector(REALSXP, ncrits));
  crits =  REAL(rcrits);
  ritercrits =  PROTECT(Rf_allocVector(REALSXP, ncrits));
  itercrits =  REAL(ritercrits);
  rdf =  PROTECT(Rf_allocVector(REALSXP, ncrits));
  df =  REAL(rdf);
  /* data in C */
  y=REAL(ry);
  S=REAL(rS);
  Dmdemi=REAL(rDmdemi);
  critsint=INTEGER(rcritsint);
  /* temp */
  gammatemp = (double *) R_alloc(n, sizeof(double));
  gammacur = (double *) R_alloc(n, sizeof(double));
  ImSkold = (double *) R_alloc(n*n, sizeof(double));
  ImSk = (double *) R_alloc(n*n, sizeof(double));
  ImS = (double *) R_alloc(n*n, sizeof(double));
  lenImSk = sizeof( double ) * n * n;
  critsold = (double *) R_alloc(ncrits, sizeof(double));
  critsok = (int *) R_alloc(ncrits, sizeof(int));
  /* init  */
  for(j = 0; j < n; j++) {
    for(i = 0; i < n; i++) {
      ImS[(j*n) + i]= -S[(j*n) + i];
    }
  }
  for(i = 0; i < n; i++) {
    gammacur[i]=y[i];
    ImS[(i*n) + i]= 1.0 + ImS[(i*n) + i];
  }
  memcpy(ImSkold, ImS, lenImSk);
  /* crit for first iteration */
  for(j = 0; j < ncrits; j++) {
    critsok[j]=0;   
    critsold[j]=MAXI;
 }
  /* loop */
  for(k = 1; k < niter; k++) {
    scr=0.0;
    for(j = 0; j < n; j++) {
      gammatemp[j] = 0.0;
    }
    /* produit: ImSk = ImSkold (I-S) */
    /* to evaluate df = tr( I - (I-S)^k ) */
    F77_NAME(dsymm)(side, uplo, &n, &n, &one, ImS, &n,
                    ImSkold, &n, &zero, ImSk, &n FCONE FCONE);
    for(i = 0; i < n; i++) {
      for(j = i; j < n; j++) {
        if (j==i) {
          gammatemp[i] = gammatemp[i] + S[(j*n) + i]*gammacur[j];
        } else {
          gammatemp[i] = gammatemp[i] + S[(j*n) + i]*gammacur[j];
          gammatemp[j] = gammatemp[j] + S[(j*n) + i]*gammacur[i];
        }      
      }
      /* SCR used in crits */
      if (symmetric) {
        scr=scr + pow(y[i] - gammatemp[i],2);
      }
      else {
        scr=scr + pow((y[i] - gammatemp[i])*Dmdemi[i],2);
      }
    }
    /* df = tr( I - (I-S)^k ) */
    /* as residuals are evaluated for iteration "k-1" */
    /* we evaluate df using ImSkold */
    dftemp=0.0;
    for (i = 0; i < n; i++) {
      dftemp=dftemp+ImSkold[(i*n) + i];
    }
    dfSold=dfS;
    dfS=n-dftemp;
    /* criteria */
    allok=1;
    for(j = 0; j < ncrits; j++) {
      if (critsok[j]==0) {
        /* crit j not already ok */
        /* calculus */
        crits[j]=calculcrit(scr, dfS, n, critsint[j]);
        if (crits[j]>critsold[j])
          {
            /* increase: crit j is ok */
            critsok[j]=1;
            /* reset to previous crit */
            crits[j] = critsold[j];
            df[j]=dfSold;
            itercrits[j]=k-1;
            /* gamma */
            for (i = 0; i < n; i++) {
              gamma[(j * n) + i]=gammacur[i];
            }
          }
        else {
          /* decrease */
          critsold[j]=crits[j];
          allok=0;
        }
      }
    }
    /* update ImSkold */
    memcpy(ImSkold, ImSk, lenImSk);
    /* test if all have already increased */
    if (allok==1) {
      break;
    } else {
      /* update gamma */
      /* gamma^{(k)} = (gamma^{(k-1)} + y) - S gamma^{(k-1)} */
      for (i = 0; i < n; i++) {
            gammacur[i]=gammacur[i] + y[i] - gammatemp[i];
      }
    }
    void R_CheckUserInterrupt(void);
  }
  /* result */
  SET_VECTOR_ELT(rans, 0, rgamma);
  SET_VECTOR_ELT(rans, 1, rdf);
  SET_VECTOR_ELT(rans, 2, rcrits);
  SET_VECTOR_ELT(rans, 3, ritercrits);
  UNPROTECT(5);
  return rans;
}

/* only product no criteria */
SEXP product(SEXP rS, /* Smoother matrix symmetric */
	         SEXP ry,
             SEXP rparamsint) 
{
  int k, i, j, n, niter, ione=1, lenImSk;
  double *S, *y, *gamma, *gammatemp, minusone=-1.0, one=1.0, zero=0.0;
  double dftemp, dfS, *ImS, *ImSkold, *ImSk;
  char *uplo="u", *side="r";
  SEXP rans, rgamma, rdf;
  /* rparamsint
     - niter number of iterations greater than 1
   */
  /* ---------- output ------------  */
  /* rans: list  */
  /*   - rgamma: vector of double (size n) D^{1/2}\beta
       - rdf: vector of double (size 1), df at niter
   */
  niter = INTEGER(rparamsint)[0];
  n=Rf_length(ry);
  /* result */
  rans=PROTECT(Rf_allocVector(VECSXP, 2));
  rgamma =  PROTECT(Rf_allocVector(REALSXP, n));
  gamma =  REAL(rgamma);
  /* data in C */
  y=REAL(ry);
  S=REAL(rS);
  /* temp */
  gammatemp = (double *) R_alloc(n, sizeof(double));
  ImSkold = (double *) R_alloc(n*n, sizeof(double));
  ImSk = (double *) R_alloc(n*n, sizeof(double));
  ImS = (double *) R_alloc(n*n, sizeof(double));
  lenImSk = sizeof( double ) * n * n;
  /* init  */
  for(j = 0; j < n; j++) {
    for(i = 0; i < n; i++) {
      ImS[(j*n) + i]= -S[(j*n) + i];
    }
  }
  for(i = 0; i < n; i++) {
    gamma[i]=y[i];
    ImS[(i*n) + i]= 1.0 + ImS[(i*n) + i];
  }
  memcpy(ImSkold, ImS, lenImSk);
  /* loop */
  if (niter>1) {
  for(k = 1; k < niter; k++) {
    for(i = 0; i < n; i++) {
      gammatemp[i] = gamma[i];
      gamma[i] = gamma[i] + y[i];
     }
    /* produit: ImSk = ImSkold (I-S) */
    F77_NAME(dsymm)(side, uplo, &n, &n, &one, ImS, &n,
                    ImSkold, &n, &zero, ImSk, &n FCONE FCONE);
    /* df = tr( I - (I-S)^k ) */
    dftemp=0.0;
    for (i = 0; i < n; i++) {
      dftemp=dftemp+ImSk[(i*n) + i];
    }
    dfS=n-dftemp;
     /* gamma^{(k)} = (gamma^{(k-1)} + y) - S gamma^{(k-1)} */
    /* result in gamma */
    F77_NAME(dsymv)(uplo, &n, &minusone, S, &n, gammatemp, &ione, &one, gamma,
                    &ione FCONE);
    /* update ImSkold */
    memcpy(ImSkold, ImSk, lenImSk);
    void R_CheckUserInterrupt(void);
  }
  } else {
    /* one iteration only ??? */
    dftemp=0.0;
    for (i = 0; i < n; i++) {
      dftemp=dftemp+ImS[(i*n) + i];
    }
    dfS=n-dftemp;
  }
  /* result */
  rdf = ScalarReal(dfS);
  SET_VECTOR_ELT(rans, 0, rgamma);
  SET_VECTOR_ELT(rans, 1, rdf);
  UNPROTECT(2);
  return rans;
}
