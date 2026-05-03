#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "kernel.h"
SEXP npreg(SEXP rx, SEXP ry, SEXP rvalx, SEXP rbw, SEXP rparamsint)
{
  int i, j;
  int nx = Rf_length(rx),  nvalx= Rf_length(rvalx), typekernel=INTEGER(rparamsint)[0];
  SEXP rans = PROTECT(Rf_allocVector(VECSXP, 2)),rregx = PROTECT(Rf_allocVector(REALSXP, nvalx));
  double some, w, wii;
  double df=0.0, bw=REAL(rbw)[0], *x=REAL(rx), *y=REAL(ry), *valx=REAL(rvalx), *regx=REAL(rregx);

  for(i = 0; i < nvalx; i++)
    regx[i] = 0.0;
  for(i = 0; i < nvalx; i++) {
	wii=0.0;
	some = 0.0;
    /* pour la i eme valeur de la grille valx :*/
    /* boucle sur les valeurs observees (indice j)*/
	for(j = 0; j < nx; j++) {
      /* poids */
      w = poidskernel(valx[i], x[j], bw, typekernel);
	  if (i==j) wii=w;
	  some=some+w;
      /* regression */ 
	  regx[i]=regx[i]+w*y[j];
	}
	if (some>0) {
	  regx[i]=regx[i]/some;
	  df=df+wii/some;
	}
  }
  /* result */
  SET_VECTOR_ELT(rans, 0, rregx);
  SET_VECTOR_ELT(rans, 1,  Rf_ScalarReal(df));
  UNPROTECT(2);
  return rans;
}

/***********************************************************/
SEXP npregcv(SEXP rx, SEXP ry, SEXP rbw, SEXP reffold, SEXP rparamsint)
{
  int i, j, k, h;
  int nx = Rf_length(rx), nbw= Rf_length(rbw), neffold= Rf_length(reffold), typekernel=INTEGER(rparamsint)[0];
  int *effold=INTEGER(reffold);
  SEXP rans = PROTECT(Rf_allocVector(VECSXP, 2)), rsse = PROTECT(Rf_allocVector(REALSXP, nbw)), rsap = PROTECT(Rf_allocVector(REALSXP, nbw));
  double *x=REAL(rx), *y=REAL(ry), *bw=REAL(rbw), *sse=REAL(rsse), *sap=REAL(rsap);
  double some, w, regx;
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
	    regx = 0.0;
	    some = 0.0;
	    /* pour la i eme valeur de la valeur en test :*/
	    /* boucle sur les valeurs en apprentissage (indice j)*/
	    for(j = 0; j < nx; j++ ) {
	      /* si pas en test */
	      if ((j>=effold[k+1])||(j<effold[k])) {
	        /* alors calculs */
	        /* poids */
            w = poidskernel(x[i], x[j], bw[h], typekernel);
	        some=some+w;
	        /* regression */ 
	        regx=regx+w*y[j];
	      }
	    } /* fin de boucle sur apprentissage*/
	    if (some>0) {
	      /* la somme par ligne */
	      regx=regx/some; 
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
