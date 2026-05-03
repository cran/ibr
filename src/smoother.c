#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h> // interrupt
#include "kernel.h"
#include "smoother.h"
SEXP Kmatrix(SEXP rx, /* Data */
             SEXP rvalx, /* values for pred */
             SEXP rbandwidth, /* bandwidths */
             SEXP rparamsint
             )
{
  int j, i, k;
  int nx, px, nvalx, typekernel,symmetric, *paramsint; 
  double *K, *bandwidth, *x, *valx;
  SEXP rK;
  x=REAL(rx);
  valx=REAL(rvalx);
  bandwidth = REAL(rbandwidth);
  /* rparamsint
     - nx number of observations
     - px number of variables
     - nvalx number of predictions
     - typekernel (1=gaussian, 2=epanechnikov, 3=quadratic, 4=uniform)
     - X and valx the same ? 1=yes and rvalx/nvalx is not used; 0=no
   */
  paramsint = INTEGER(rparamsint);
  nx=paramsint[0];
  px=paramsint[1];
  nvalx=paramsint[2];
  typekernel=paramsint[3];
  symmetric=paramsint[4];
  /* ---------- output ------------  */
  /* rK: matrix of double K (nvalx x nx)
   */
  rK = PROTECT(Rf_allocMatrix(REALSXP, nvalx, nx));
  K=REAL(rK);
  /* case symmetric x=valx */
  if (symmetric==1) {  
    for(i = 0; i < nx; i++) {
      for (j= i; j < nx; j++) { 
        K[(nx*j)+i]= 1.0;
        for (k=0; k< px;k++) {
	      K[(nx*j)+i]= K[(nx*j)+i]*poidskernel(x[k*nx+i], x[k*nx+j],bandwidth[k],typekernel);
        }
        K[(nx*i)+j]=K[(nx*j)+i];
      }
      void R_CheckUserInterrupt(void);
    }
  } else {
    /* case non symmetric x .ne. valx */
    for(i = 0; i < nvalx; i++) {
	  for (j= 0; j < nx; j++) { 
        K[(nvalx*j)+i]= 1.0;
        for (k=0; k< px;k++) {
	      K[(nvalx*j)+i]= K[(nvalx*j)+i]*poidskernel(valx[k*nvalx+i], x[k*nx+j],bandwidth[k],typekernel);
        }
      }
      void R_CheckUserInterrupt(void); 
    }
  }
  UNPROTECT(1);
  return rK;
}

SEXP Smatrix(SEXP rK) /* K matrix */
{
  int j, i, nK, pK;
  double *K, rowsum, *S;
  SEXP rS;
  pK=Rf_ncols(rK);
  nK=Rf_nrows(rK);
  K=REAL(rK);
  /* ---------- output ------------  */
  /* rS: matrix of double K (nvalx x nx)
   */
  rS = PROTECT(Rf_allocMatrix(REALSXP, nK, pK));
  S=REAL(rS);
  for(i = 0; i < nK; i++) {
      rowsum=0.0;
      for (j= 0; j < pK; j++) {
        rowsum=rowsum + K[(nK*j)+i];
      }
      for (j= 0; j < pK; j++) {
        S[(nK*j)+i]=K[(nK*j)+i]/rowsum;
      }
      void R_CheckUserInterrupt(void);
    }
   UNPROTECT(1);
  return rS;
}

SEXP Amatrix(SEXP rK) /* K matrix */
{
  int j, i, nK, pK;
  double *K, *rowsums, *colsums, *A;
  SEXP rans, rA, rrowsums, rcolsums;
  pK=Rf_ncols(rK);
  nK=Rf_nrows(rK);
  K=REAL(rK);
  /* ---------- output ------------  */
  /* rans: list  */
  /*   - rS: matrix of double D_r^{-1/2} K D_c^{-1/2}
       - rrowsums = diag(D_r) vector of length nK
       - colsums = diag(D_c) vector of length pK
   */
  rans=PROTECT(Rf_allocVector(VECSXP, 3));
  rA = PROTECT(Rf_allocMatrix(REALSXP, nK, pK));
  rrowsums = PROTECT(Rf_allocVector(REALSXP, nK));
  rcolsums = PROTECT(Rf_allocVector(REALSXP, pK));
  A=REAL(rA);
  rowsums = REAL(rrowsums);
  colsums = REAL(rcolsums);
  
  /* rowsums and colsums*/
  for (j= 0; j < pK; j++) {
      colsums[j]=0.0;
  }
  for(i = 0; i < nK; i++) {
    rowsums[i]=0.0;
    for (j= 0; j < pK; j++) {
      rowsums[i]=rowsums[i] + K[(nK*j)+i];
      colsums[j]=colsums[j] + K[(nK*j)+i];
    }
    void R_CheckUserInterrupt(void);
  }
  for(i = 0; i < nK; i++) {
    for (j= 0; j < pK; j++) {
      A[(nK*j)+i] = K[(nK*j)+i]/sqrt(rowsums[i])/sqrt(colsums[j]);
    }
  }
   /* result */
  SET_VECTOR_ELT(rans, 0, rA);
  SET_VECTOR_ELT(rans, 1, rrowsums);
  SET_VECTOR_ELT(rans, 2, rcolsums);
  UNPROTECT(4);
  return rans;
}

