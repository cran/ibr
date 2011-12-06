#include <float.h> /* DBL_EPSILON */
#include <math.h>
#include <R.h>
#include <R_ext/Linpack.h>	/* Fortran routines */

/* Copyright (C) pac@univ-rennes2.fr

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License   
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
(www.gnu.org/copyleft/gpl.html)

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
USA.*/

void semikerlog(double *x, double *xk, int *nx, int *nxk,double *ksurdeux, int *px, int *negatif, double *dista, int *symmetric)
{
  double eps;
  int i,j,l,k;
  eps = DBL_EPSILON;
  /* case symmetric x=xk */
  if (*symmetric) {
  for(i = 0; i < *nx; i++) {
    dista[((*nx)*i)+i]= 0.0;
    for (j= i; j < *nx; j++) { 
      /* computing squared euclidean distance */
      dista[((*nx)*j)+i]= 0.0;
      for (l=0; l< *px;l++) {
 	dista[((*nx)*j)+i]= dista[((*nx)*j)+i]+
	  pow( (x[l*(*nx)+i]-x[l*(*nx)+j]),2);
      }
      /* semi kernel : distance^k * log(distance)*/
      if ( dista[((*nx)*j)+i]< eps) {
	dista[((*nx)*j)+i]= 0.0; 
      } else {
	/* changement de signe pour etre coherent avec le noyau Tps
	   dont on connait la constante */
	if (*negatif) {
	  dista[((*nx)*j)+i]= - pow(dista[((*nx)*j)+i],*ksurdeux) * log(dista[((*nx)*j)+i])/2; 
	} else {
	  dista[((*nx)*j)+i]= pow(dista[((*nx)*j)+i],*ksurdeux) * log(dista[((*nx)*j)+i])/2; 
	}
      }
	/* symmetric */
      dista[((*nx)*i)+j]=dista[((*nx)*j)+i];
    }
  }
  } else {
  /* case non symmetric x .ne. xk */
  for(i = 0; i < *nx; i++) {
    for (j= 0; j < *nxk; j++) { 
      /* computing squared euclidean distance */
      dista[((*nx)*j)+i]= 0.0;
      for (l=0; l< *px;l++) {
 	dista[((*nx)*j)+i]= dista[((*nx)*j)+i]+
	  pow( (x[l*(*nx)+i]-xk[l*(*nxk)+j]),2);
      }
      /* semi kernel : distance^k * log(distance)*/
      if ( dista[((*nx)*j)+i]< eps) {
	dista[((*nx)*j)+i]= 0.0; 
      } else {
	/* changement de signe pour etre coherent avec le noyau Tps
	   dont on connait la constante */
	if (*negatif) {
	  dista[((*nx)*j)+i]= - pow(dista[((*nx)*j)+i],*ksurdeux) * log(dista[((*nx)*j)+i])/2; 
	} else {
	  dista[((*nx)*j)+i]= pow(dista[((*nx)*j)+i],*ksurdeux) * log(dista[((*nx)*j)+i])/2; 
	}
      }
    }
  }
  }
}
void semikerpow(double *x, double *xk, int *nx, int *nxk,double *ksurdeux, int *px, int *negatif, double *dista, int *symmetric)
{
  double eps;
  int i,j,l,k;
  eps = DBL_EPSILON;
  /* case symmetric x=xk */
  if (*symmetric) {
  for(i = 0; i < *nx; i++) {
    dista[((*nx)*i)+i]= 0.0;
    for (j= i; j < *nx; j++) { 
      /* computing squared euclidean distance */
      dista[((*nx)*j)+i]= 0.0;
      for (l=0; l< *px;l++) {
 	dista[((*nx)*j)+i]= dista[((*nx)*j)+i]+
	  pow( (x[l*(*nx)+i]-x[l*(*nx)+j]),2);
      }
      /* semi kernel : distance^k * log(distance)*/
      if ( dista[((*nx)*j)+i] < eps) {
	dista[((*nx)*j)+i]= 0.0; 
      } else {
	/* changement de signe pour etre coherent avec le noyau Tps
	   dont on connait la constante */
	if (*negatif) {
	  dista[((*nx)*j)+i]= - pow(dista[((*nx)*j)+i],*ksurdeux) ;
	} else {
	  dista[((*nx)*j)+i]= pow(dista[((*nx)*j)+i],*ksurdeux) ;
	}
      }
	/* symmetric */
      dista[((*nx)*i)+j]=dista[((*nx)*j)+i];
    }
  }
  } else {
  /* case non symmetric x .ne. xk */
  for(i = 0; i < *nx; i++) {
    for (j= 0; j < *nxk; j++) { 
      /* computing squared euclidean distance */
      dista[((*nx)*j)+i]= 0.0;
      for (l=0; l< *px;l++) {
 	dista[((*nx)*j)+i]= dista[((*nx)*j)+i]+
	  pow( (x[l*(*nx)+i]-xk[l*(*nxk)+j]),2);
      }
      /* semi kernel : distance^k * log(distance)*/
      if ( dista[((*nx)*j)+i] < eps) {
	dista[((*nx)*j)+i]= 0.0; 
      } else {
	/* changement de signe pour etre coherent avec le noyau Tps
	   dont on connait la constante */
	if (*negatif) {
	  dista[((*nx)*j)+i]= - pow(dista[((*nx)*j)+i],*ksurdeux) ;
	} else {
	  dista[((*nx)*j)+i]= pow(dista[((*nx)*j)+i],*ksurdeux) ;
	}
      }
    }
  }
  }
}
void polynom(int *m, int *n, int *dim, double 
	*des, int *lddes, int *npoly, double *t, int *ldt, 
	   int *wptr, int *info)
{
    /* System generated locals */
  int one;
  one=1;
  
  /* Local variables */
  static int i, j, k, jj, ii, nt, tt, bptr, eptr;
  
  
  /* Purpose: create t matrix and append s1 to it. */
  
  /* On Entry: */
  /*   m			order of the derivatives in the penalty */
  /*   n			number of rows in des */
  /*   dim			number of columns in des */
  /*   des(lddes,dim)	variables to be splined */
  /*   lddes		leading dimension of des as declared in the */
  /* 			calling program */
  /*   ldt			leading dimension of t as declared in the */
  /* 			calling program */
  
  /*   npoly		dimension of polynomial part of spline */
  /* On Exit: */
  /*   t(ldt,npoly+ncov1)	[t:s1] */
  /*   info 		error indication */
  /*   			   0 : successful completion */
  /* 		 	   1 : error in creation of t */
  /* Work Arrays: */
  /*   wptr(dim)		integer work vector */
  
  /* Subprograms Called Directly: */
  /* 	Blas  - dcopy */
  /* 	Other - No */
  /* Function Body */
  *info = 0;
  /*      npoly = mkpoly(m,dim) */
  /* constant */
  for (j = 0; j < *n; ++j) {
    t[j] = 1.0;
  }
  /* init nt */
  nt = 0;
  if (*npoly > 1) {
  /*  X*/
    for (j = 0; j < *dim; ++j) {
      nt = j + 1;
      wptr[j] = nt;
      /* copy x */
       F77_CALL(dcopy)(n, &des[j * *lddes ], &one, &t[nt * *ldt], &one);
    }
    
    /*     get cross products of x's in null space for m>2 */
    
    /*     WARNING: do NOT change next do loop unless you fully understand: */
    /*              This first gets x1*x1, x1*x2, x1*x3, then */
    /*              x2*x2, x2*x3, and finally x3*x3 for dim=3,n=3 */
    /*              wptr(1) is always at the beginning of the current */
    /* 	       level of cross products, hence the end of the */
    /* 	       previous level which is used for the next. */
    /* 	       wptr(j) is at the start of xj * (previous level) */
    /* loop on degree */
    if (*m > 2) {
    for (k = 2; k <= (*m - 1); ++k) {
      /* loop on variable */
      for (j = 0; j < *dim; ++j) {
	/* no need to make all products */
	bptr = wptr[j];
	wptr[j] = nt + 1;
	eptr = wptr[0] - 1;
	for (tt = bptr; tt <= eptr; ++tt) {
	  ++nt;
	  for (ii = 0; ii < *n; ++ii) {
	    t[ii + nt * *ldt] = des[ii + j * *lddes] * t[ii + tt * *ldt];
	  }
	}
      }
    }
    }
    if (nt != (*npoly-1)) {
      *info = 1;
    }
  }
} 
