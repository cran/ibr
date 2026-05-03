#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "kernel.h"

/* Regression noyau gaussien */
SEXP npregpol(SEXP rx, SEXP ry, SEXP rvalx, SEXP rbw, SEXP rparamsint)
{
  int i, j ; 
  int nx = Rf_length(rx),  nvalx= Rf_length(rvalx), typekernel=INTEGER(rparamsint)[0];
  SEXP rans = PROTECT(Rf_allocVector(VECSXP, 3)), rregx = PROTECT(Rf_allocVector(REALSXP, nvalx)), rderiv = PROTECT(Rf_allocVector(REALSXP, nvalx));
  double S0, S1, S2, w, T0, T1, wii;
  double df, bw=REAL(rbw)[0], *x=REAL(rx), *y=REAL(ry), *valx=REAL(rvalx), *regx=REAL(rregx), *deriv=REAL(rderiv);
  /* initialisation */      
  w = 0.0;
  df=0.0;
  for(i = 0; i < nvalx; i++)
    regx[i] = 0.0;
  for(i = 0; i < nvalx; i++) {
    wii=0.0;
    S0 = 0.0;
    S1 = 0.0;
    S2 = 0.0;
    T1 = 0.0;
    T0 = 0.0;
    /* pour la i eme valeur de la grille valx :*/
    /* boucle sur les valeurs observees (indice j)*/
    for(j = 0; j < nx; j++) {
      /* poids */
      w = poidskernel(valx[i], x[j], bw, typekernel);
      if (i==j) wii=w;
      S0 = S0 + w;
      S1 = S1 + w * (x[j] - valx[i]);
      S2 = S2 + w * pow((x[j] - valx[i]),2);
      /* regression */ 
      T0 = T0+w*y[j];
      T1 = T1+ (x[j] - valx[i]) * w * y[j];
    }
    if (S0>0) {
      regx[i]= (S2 * T0 - S1 * T1)/(S0 * S2 - pow(S1,2));
      deriv[i]=(- S1 * T0 + S0 * T1)/(S0 * S2 - pow(S1,2));
      df=df+wii/S0;
    }
  }
  /* result */
  SET_VECTOR_ELT(rans, 0, rregx);
  SET_VECTOR_ELT(rans, 1, rderiv);
  SET_VECTOR_ELT(rans, 2,  Rf_ScalarReal(df));
  UNPROTECT(3);
  return rans;
}

/***********************************************************/
SEXP npregpolcv(SEXP rx, SEXP ry, SEXP rbw, SEXP reffold, SEXP rparamsint)
{
  int i, j , k, h ; 
  int nx = Rf_length(rx), nbw= Rf_length(rbw), neffold= Rf_length(reffold), typekernel=INTEGER(rparamsint)[0];
  int *effold=INTEGER(reffold);
  SEXP rans = PROTECT(Rf_allocVector(VECSXP, 2)), rsse = PROTECT(Rf_allocVector(REALSXP, nbw)), rsap = PROTECT(Rf_allocVector(REALSXP, nbw));
  double *x=REAL(rx), *y=REAL(ry), *bw=REAL(rbw), *sse=REAL(rsse), *sap=REAL(rsap);
  double S0, S1, S2, w, T0, T1,  regx;
  /* initialisation */      
  w = 0.0;
  neffold=neffold-1;
  /* boucle sur les fenetres*/
  for (h=0; h < nbw; h++) {
    sse[h]=0.0;
    sap[h]=0.0;
    /* boucle sur les fold */
    for (k=0; k < neffold; k++) {
      /* boucle sur la partie test */
      for(i = effold[k]; i < effold[k+1]; i++) {
	    S0 = 0.0;
	    S1 = 0.0;
	    S2 = 0.0;
	    T1 = 0.0;
	    T0 = 0.0;
	    /* pour la i eme valeur de la grille x :*/
	    /* boucle sur les valeurs observees (indice j)*/
	    for(j = 0; j < nx; j++) {
	      /* si pas en test */
	      if ((j>=effold[k+1])||(j<effold[k])) {
	        /* poids */
            w = poidskernel(x[i], x[j], bw[h], typekernel);
	        S0 = S0 + w;
	        S1 = S1 + w * (x[j] - x[i]);
	        S2 = S2 + w * pow((x[j] - x[i]),2);
	        /* regression */ 
	        T0 = T0+w*y[j];
	        T1 = T1+ (x[j] - x[i]) * w * y[j];
	      }
	    } /* fin de boucle sur apprentissage*/
	    if (S0>0) {
	      regx= (S2 * T0 - S1 * T1)/(S0 * S2 - pow(S1,2));
	      /*les ecarts pour la fenetre h*/
	      sse[h]=sse[h]+pow(y[i]-regx,2);
	      sap[h]=sap[h]+fabs((y[i]-regx)/y[i]);
	    } else {
	      sse[h]=sse[h]+pow(y[i],2);
	      sap[h]=sap[h]+1;
	    }
      }
    }
  }
  /* result */
  SET_VECTOR_ELT(rans, 0, rsse);
  SET_VECTOR_ELT(rans, 1, rsap);
  UNPROTECT(3);
  return rans;
}
