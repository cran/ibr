/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999, 2001 the R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

/* from NETLIB c/brent.shar with max.iter, add'l info and convergence
   details hacked in by Peter Dalgaard */
/* modified for calculating bandwidths for given effective df*/

/*************************************************************************
 *			    C math library
 * function ZEROIN - obtain a function zero within the given range
 *
 * Input
 *	double zeroin(ax,bx,f,info,Tol,Maxit)
 *	double ax;			Root will be seeked for within
 *	double bx;			a range [ax,bx]
 *	double (*f)(double x, void *info); Name of the function whose zero
 *					will be seeked for
 *	void *info;			Add'l info passed to f
 *	double *Tol;			Acceptable tolerance for the root
 *					value.
 *					May be specified as 0.0 to cause
 *					the program to find the root as
 *					accurate as possible
 *
 *	int *Maxit;			Max. iterations
 *
 *
 * Output
 *	Zeroin returns an estimate for the root with accuracy
 *	4*EPSILON*abs(x) + tol
 *	*Tol returns estimated precision
 *	*Maxit returns actual # of iterations, or -1 if maxit was
 *      reached without convergence.
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.180 of the Russian edition
 *
 *	The function makes use of the bisection procedure combined with
 *	the linear or quadric inverse interpolation.
 *	At every step program operates on three abscissae - a, b, and c.
 *	b - the last and the best approximation to the root
 *	a - the last but one approximation
 *	c - the last but one or even earlier approximation than a that
 *		1) |f(b)| <= |f(c)|
 *		2) f(b) and f(c) have opposite signs, i.e. b and c confine
 *		   the root
 *	At every step Zeroin selects one of the two new approximations, the
 *	former being obtained by the bisection procedure and the latter
 *	resulting in the interpolation (if a,b, and c are all different
 *	the quadric interpolation is utilized, otherwise the linear one).
 *	If the latter (i.e. obtained by the interpolation) point is
 *	reasonable (i.e. lies within the current interval [b,c] not being
 *	too close to the boundaries) it is accepted. The bisection result
 *	is used in the other case. Therefore, the range of uncertainty is
 *	ensured to be reduced at least by the factor 1.6
 *
 ************************************************************************
 */
#include <R.h>
#include <Rinternals.h>
#include "kernel.h"
static double * vector(int nch)
{
    double *m;
    m = (double*) R_alloc((nch + 1), sizeof(double));
    return m;
}



/* Trace noyau gaussien */
static double caltrkernel(double bw, double *x, int *nx, double *objectif, int colonne, double *sl, int typekernel) 
{
  int j, i; 
  double  tmp, trtmp, trace; 
  trace=0.0;
  trtmp=0.0;
  for(i = 0; i < *nx; i++) {
	sl[i]=0.0;
  }
  for(i = 0; i < *nx; i++) {
    for (j= i; j < *nx; j++) {
      tmp=poidskernel(x[colonne*(*nx)+i], x[colonne*(*nx)+j], bw, typekernel);
      sl[i]=sl[i]+tmp;
      if (j==i)
	trtmp= tmp;
      else 
	sl[j]=tmp+sl[j];
    }
    trace=trace+trtmp/sl[i];
  }
  trace=trace- *objectif;
  return trace;
}


#define EPSILON DBL_EPSILON

void zerotrace(			/* An estimate of the root */
               double *ax,				/* Left border | of the range	*/
               double *bx,			/* Right border| the root is seeked*/
               /*    double *fa, double *fb,		 f(a), f(b) [now calculated] */
               double *x, int *nx, int *px, double *objectif,    /* Arguments for trace*/
               double *Tol,			/* Acceptable tolerance		*/
               int *Maxit, double *bandwidth, /* Max # of iterations */
               int typekernel)				
{
  double a, b, c, fc, fa, fb;			/* Abscissae, descr. see above,  f(c) */
  double tol, *sl;
  int maxit, k, success,pb;
    sl=vector(*nx);
    pb=0;
    for (k = 0; k < *px; k++) {
      success=0;
      b=*bx;
      a = ax[k];  
      fa = caltrkernel(a,x,nx,objectif,k,sl,typekernel);
      if(fa == 0.0) {
	bandwidth[k]=a;
	/*	return a;*/
	success=1;
	break;
      }
      while(fa>0) {
	a = a * 2;
	fa = caltrkernel(a,x,nx,objectif,k,sl,typekernel);
      }
      if(fa == 0.0) {
	bandwidth[k]=a;
	/*	return a;*/
	success=1;
	break;
      }
      
      fb = caltrkernel(b,x,nx,objectif,k,sl,typekernel);
      if(fb ==  0.0) {
	bandwidth[k]=b;
	success=1;
	break;
      }
      while(fb<0) {
	b = b / 2;
	fb = caltrkernel(b,x,nx,objectif,k,sl,typekernel);
      }
      if(fb ==  0.0) {
	bandwidth[k]=b;
	success=1;
	break;
      }
      c = a;   fc = fa;
      maxit = *Maxit + 1; tol = *Tol;
      /* First test if we have found a root at an endpoint */

    while(maxit--)		/* Main iteration loop	*/
    {
	double prev_step = b-a;		/* Distance from the last but one
					   to the last approximation	*/
	double tol_act;			/* Actual tolerance		*/
	double p;			/* Interpolation step is calcu- */
	double q;			/* lated in the form p/q; divi-
					 * sion operations is delayed
					 * until the last moment	*/
	double new_step;		/* Step at this iteration	*/

	if( fabs(fc) < fabs(fb) )
	{				/* Swap data for b to be the	*/
	    a = b;  b = c;  c = a;	/* best approximation		*/
	    fa=fb;  fb=fc;  fc=fa;
	}
	tol_act = 2*EPSILON*fabs(b) + tol/2;
	new_step = (c-b)/2;

	if( fabs(new_step) <= tol_act || fb == (double)0 )
	{
	    bandwidth[k]=b;
	    success=1;
	    break ;			/* Acceptable approx. is found	*/
	}

	/* Decide if the interpolation can be tried	*/
	if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
	    && fabs(fa) > fabs(fb) ) {	/* and was in true direction,
					 * Interpolation may be tried	*/
	    register double t1,cb,t2;
	    cb = c-b;
	    if( a==c ) {		/* If we have only two distinct	*/
					/* points linear interpolation	*/
		t1 = fb/fa;		/* can only be applied		*/
		p = cb*t1;
		q = 1.0 - t1;
	    }
	    else {			/* Quadric inverse interpolation*/

		q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
		p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
		q = (q-1.0) * (t1-1.0) * (t2-1.0);
	    }
	    if( p>(double)0 )		/* p was calculated with the */
		q = -q;			/* opposite sign; make p positive */
	    else			/* and assign possible minus to	*/
		p = -p;			/* q				*/

	    if( p < (0.75*cb*q-fabs(tol_act*q)/2) /* If b+p/q falls in [b,c]*/
		&& p < fabs(prev_step*q/2) )	/* and isn't too large	*/
		new_step = p/q;			/* it is accepted
						 * If p/q is too large then the
						 * bisection procedure can
						 * reduce [b,c] range to more
						 * extent */
	}

	if( fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
	    if( new_step > (double)0 )	/* than tolerance		*/
		new_step = tol_act;
	    else
		new_step = -tol_act;
	}
	a = b;	fa = fb;			/* Save the previous approx. */
	b += new_step;	fb = caltrkernel(b,x,nx,objectif,k,sl,typekernel);	/* Do step to a new approxim. */
	if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
	    /* Adjust c for it to have a sign opposite to that of b */
	    c = a;  fc = fa;
	}

    }
    /* failed! */
    if (success==0) {
      bandwidth[k]=b;
      pb = 1;
    }
    }
    if (pb==1) *Maxit=-1;
}


void evaltracetotal(double *x, int *nx, int *px,
			     double *ax, double *bx, double *objectifuni,
		     double *Tol, int *Maxit, double *objectif, 
		            double *resubw, double *trtot, double *sl, int typekernel) /*results */	
{
  double a, b, c, fc, fa, fb;			/* Abscissae, descr. see above,  f(c) */
  double tol;
    int maxit, k, success, i ,j;
    for (k = 0; k < *px; k++) {
      b=bx[k];
      success=0;
      a = ax[k];
      fb = caltrkernel(b,x,nx,objectifuni,k,sl,typekernel);
      fa = caltrkernel(a,x,nx,objectifuni,k,sl,typekernel);
      maxit = *Maxit+ 1; tol = *Tol;
      if(fa == 0.0) {
	resubw[k]=a;
	/*	return a;*/
	success=1;
	maxit=0;
    }
      while(fa>0) {
	a = a * 2;
	fa = caltrkernel(a,x,nx,objectifuni,k,sl,typekernel);
      }

    if(fb ==  0.0) {
	resubw[k]=b;
	success=1;
	maxit=0;
    }

      while(fb<0) {
	b = b / 2;
	fb = caltrkernel(b,x,nx,objectifuni,k,sl,typekernel);
      }
      /* depart */
      c = a;   fc = fa;
    while(maxit--)		/* Main iteration loop	*/
    {
	double prev_step = b-a;		/* Distance from the last but one
					   to the last approximation	*/
	double tol_act;			/* Actual tolerance		*/
	double p;			/* Interpolation step is calcu- */
	double q;			/* lated in the form p/q; divi-
					 * sion operations is delayed
					 * until the last moment	*/
	double new_step;		/* Step at this iteration	*/

	if( fabs(fc) < fabs(fb) )
	{				/* Swap data for b to be the	*/
	    a = b;  b = c;  c = a;	/* best approximation		*/
	    fa=fb;  fb=fc;  fc=fa;
	}
	tol_act = 2*EPSILON*fabs(b) + tol/2;
	new_step = (c-b)/2;

	if( fabs(new_step) <= tol_act || fb == (double)0 )
	{
	    resubw[k]=b;
	    success=1;
	    break ;			/* Acceptable approx. is found	*/
	}

	/* Decide if the interpolation can be tried	*/
	if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
	    && fabs(fa) > fabs(fb) ) {	/* and was in true direction,
					 * Interpolation may be tried	*/
	    register double t1,cb,t2;
	    cb = c-b;
	    if( a==c ) {		/* If we have only two distinct	*/
					/* points linear interpolation	*/
		t1 = fb/fa;		/* can only be applied		*/
		p = cb*t1;
		q = 1.0 - t1;
	    }
	    else {			/* Quadric inverse interpolation*/

		q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
		p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
		q = (q-1.0) * (t1-1.0) * (t2-1.0);
	    }
	    if( p>(double)0 )		/* p was calculated with the */
		q = -q;			/* opposite sign; make p positive */
	    else			/* and assign possible minus to	*/
		p = -p;			/* q				*/

	    if( p < (0.75*cb*q-fabs(tol_act*q)/2) /* If b+p/q falls in [b,c]*/
		&& p < fabs(prev_step*q/2) )	/* and isn't too large	*/
		new_step = p/q;			/* it is accepted
						 * If p/q is too large then the
						 * bisection procedure can
						 * reduce [b,c] range to more
						 * extent */
	}

	if( fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
	    if( new_step > (double)0 )	/* than tolerance		*/
		new_step = tol_act;
	    else
		new_step = -tol_act;
	}
	a = b;	fa = fb;			/* Save the previous approx. */
	b += new_step;	
	fb = caltrkernel(b,x,nx,objectifuni,k,sl, typekernel);/* Do step to a new approxim. */
	if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
	    /* Adjust c for it to have a sign opposite to that of b */
	    c = a;  fc = fa;
	}

    }
    /* failed! */
    if (success==0) {
      resubw[k]=b;
    }
    }
    /* maintenant calcul de la trace totale */    

  double  tmp, trtmp, trace; 
  trace=0.0;
  trtmp=0.0;
  for(i = 0; i < *nx; i++) {
	sl[i]=0.0;
  }
  for(i = 0; i < *nx; i++) {
    for (j= i; j < *nx; j++) {
      tmp=1.0;
      for (k=0; k< *px;k++) {
	    tmp= tmp*poidskernel(x[k*(*nx)+i], x[k*(*nx)+j], resubw[k], typekernel);
      }
      sl[i]=sl[i]+tmp;
      if (j==i)
	trtmp= tmp;
      else 
	sl[j]=tmp+sl[j];
    }
    trace=trace+trtmp/sl[i];
  }
  *trtot=trace- *objectif;
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

void zerotracetotal(double *ax, 
                    /* les fenetres etroites et resultats 1ere coord. inutile*/
                    double *bx,	/* les fenetres large 1ere coord. inutile*/
                    double *x, int *nx, int *px, /* Donnes */
			        /*    double *objectifuni,   objectif univarie*/
			        double *Tol,  /* Acceptable tolerance*/
			        int *Maxit, /* Max  of iterations*/ 
			        double *objectif, /*objectif trace totale*/
			        double *output, /*output*/
                    int typekernel)
{
  double a, b, c, fc, fa, fb;	/* Abscissae, descr. see above,  f(c) */
  double tol,tmp, *sl ,trtot, fpx;
  int maxit,k;
    sl=vector(*nx);
    tmp=0.0;
    maxit=*Maxit + 1;
    tol=*Tol;
    fpx= (double)(*px) ;
    tmp=0.0;
    k=0;
    a=pow(*objectif,1.0/fpx);
    /*evaluation de l'objectif */
    evaltracetotal(x,nx,px,ax,bx,&a,Tol,Maxit,objectif,output,&trtot,sl, typekernel);
    fa = trtot;
    if(fa ==  0.0) {
      *Maxit=0;
      return;
    }
    if(fa < 0.0) {
      /* bornes en ax ok */  
      for (k = 0; k < (*px); k++) {
	ax[k]=output[k];
      }
      b=pow(*objectif,2.0/fpx);
      evaltracetotal(x,nx,px,ax,bx,&b,Tol,Maxit,objectif,output,&trtot,sl, typekernel);
      fb = trtot;
      if(fb ==  0.0) {
	*Maxit=0;
	return;
      }
      if (fb<0.0) {
	/* grande trace (grand objectif), fenetre petite */
	b=*objectif; 
	evaltracetotal(x,nx,px,ax,bx,&b,Tol,Maxit,objectif,
		           output,&trtot,sl, typekernel);
	fb = trtot;
	if(fb ==  0.0) {
	  *Maxit=0;
	  return;
	}
	if (fb<0.0) {
	  b=*nx-2; 
	  evaltracetotal(x,nx,px,ax,bx,&b,Tol,Maxit,objectif,
			         output,&trtot,sl, typekernel);
	  fb = trtot;
	  if(fb ==  0.0) {
	    *Maxit=0;
	    return;
	  } 
	  if (fb<0.0) {
	    *Maxit=-1;
	    return;
	  }
	}
      }
      /* afiner la borne  */
      for (k = 0; k < (*px); k++) {
	bx[k]=output[k];
      }
    } 
    else {
      /* en theorie impossible !*/
      fb=fa;
      b=a;
      /* bornes en bx ok */  
      for (k = 0; k < (*px); k++) {
	bx[k]=output[k];
      }
      a=pow(*objectif,1.0/(2*fpx));
      evaltracetotal(x,nx,px,ax,bx,&a,Tol,Maxit,objectif,output,&trtot,sl, typekernel);
      fa = trtot;
      if(fa ==  0.0) {
	*Maxit=0;
	return;
      }     
      if (fa > 0) {
	k=0;
	a=caltrkernel(ax[k],x,nx,&tmp,k,sl,typekernel);
	evaltracetotal(x,nx,px,ax,bx,&a,Tol,Maxit,objectif,output,&trtot,sl,typekernel);
	fa = trtot;
	if(fa ==  0.0) {
	  *Maxit=0;
	  return;
	}     
	if (fa > 0) {
	  for (k = 0; k < (*px); k++) {
	    ax[k]=4*ax[k];
	  }
	  k=0;
	  a=caltrkernel(ax[k],x,nx,&tmp,k,sl,typekernel);
	  evaltracetotal(x,nx,px,ax,bx,&a,Tol,Maxit,objectif,output,&trtot,sl,typekernel);
	  fa = trtot;
	  if(fa ==  0.0) {
	    *Maxit=0;
	    return;
	  }
	  if (fa > 0) {
	    *Maxit=-1;
	    return;
	  }
	}
      }
      /* afiner la borne  */
      for (k = 0; k < (*px); k++) {
	ax[k]=output[k];
      }
    }
    /* mise en route */
    c = a;   fc = fa;
    
    while(maxit--)		/* Main iteration loop	*/
    {
	double prev_step = b-a;		/* Distance from the last but one
					   to the last approximation	*/
	double tol_act;			/* Actual tolerance		*/
	double p;			/* Interpolation step is calcu- */
	double q;			/* lated in the form p/q; divi-
					 * sion operations is delayed
					 * until the last moment	*/
	double new_step;		/* Step at this iteration	*/

	if( fabs(fc) < fabs(fb) )
	{				/* Swap data for b to be the	*/
	    a = b;  b = c;  c = a;	/* best approximation		*/
	    fa=fb;  fb=fc;  fc=fa;
	}
	tol_act = 2*EPSILON*fabs(b) + tol/2;
	new_step = (c-b)/2;

	if( fabs(new_step) <= tol_act || fb == (double)0 )
	{
	  *Maxit -= maxit;
	  return;  /* Acceptable approx. is found	*/
	}

	/* Decide if the interpolation can be tried	*/
	if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
	    && fabs(fa) > fabs(fb) ) {	/* and was in true direction,
					 * Interpolation may be tried	*/
	    register double t1,cb,t2;
	    cb = c-b;
	    if( a==c ) {		/* If we have only two distinct	*/
					/* points linear interpolation	*/
		t1 = fb/fa;		/* can only be applied		*/
		p = cb*t1;
		q = 1.0 - t1;
	    }
	    else {			/* Quadric inverse interpolation*/

		q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
		p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
		q = (q-1.0) * (t1-1.0) * (t2-1.0);
	    }
	    if( p>(double)0 )		/* p was calculated with the */
		q = -q;			/* opposite sign; make p positive */
	    else			/* and assign possible minus to	*/
		p = -p;			/* q				*/

	    if( p < (0.75*cb*q-fabs(tol_act*q)/2) /* If b+p/q falls in [b,c]*/
		&& p < fabs(prev_step*q/2) )	/* and isn't too large	*/
		new_step = p/q;			/* it is accepted
						 * If p/q is too large then the
						 * bisection procedure can
						 * reduce [b,c] range to more
						 * extent */
	}

	if( fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
	    if( new_step > (double)0 )	/* than tolerance		*/
		new_step = tol_act;
	    else
		new_step = -tol_act;
	}
	a = b;	fa = fb;			/* Save the previous approx. */
	b += new_step;
	evaltracetotal(x,nx,px,ax,bx,&b,Tol,Maxit,objectif,output,&trtot,sl,typekernel);
	fb = trtot;
	if (fb<0) {
	  /* sauvegarde borne etc */
	  for (k = 0; k < (*px); k++) {
	    ax[k]=output[k];
	  }
	}
	else {
	  /* sauvegarde borne etc */
	  for (k = 0; k < (*px); k++) {
	    bx[k]=output[k];
	  }
	}	  
	if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
	    /* Adjust c for it to have a sign opposite to that of b */
	    c = a;  fc = fa;
	}

    }
    /* failed! */
    *Maxit = -1;
}

/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/
/*************************************************************************/

SEXP choosebw(SEXP rax, /* les fenetres etroites et resultats 1ere coord. inutile*/
	           SEXP rbx,	/* les fenetres large 1ere coord. inutile*/
	           SEXP rx, /* Data */
               SEXP rparamsdouble,
               SEXP rparamsint
               ) 
{
  int nx, px, Maxit, dftotal, typekernel, *paramsint; 
  double *bandwidth, Tol, objectif, *paramsdouble, *ax, *bx, *x;
  SEXP rbandwidth;
  x=REAL(rx);
  ax=REAL(rax);
  bx=REAL(rbx);
  /* rparamsdouble
     - Acceptable tolerance
     - total trace objective
   */
  paramsdouble =  REAL(rparamsdouble);
  Tol=paramsdouble[0];
  objectif=paramsdouble[1];
  /* rparamsint
     - nx number of observations
     - px number of variables
     - Max  of iterations
     - dftotal (if one then objective is total df)
     - typekernel (1=gaussian, 2=epanechnikov, 3=quadratic, 4=uniform)
   */
  paramsint = INTEGER(rparamsint);
  nx=paramsint[0];
  px=paramsint[1];
  Maxit=paramsint[2];
  dftotal=paramsint[3];
  typekernel=paramsint[4];
  /* ---------- output ------------  */
  /* bandwidth (vector of double of length px: bandwidth for each col/variable)
   */
  rbandwidth = PROTECT(Rf_allocVector(REALSXP, px));
  bandwidth = REAL(rbandwidth);

  if (dftotal==1) {
    zerotracetotal(ax,bx,x,&nx,&px,&Tol,&Maxit,&objectif,bandwidth,typekernel);
  } 
  else {
    zerotrace(ax,bx,x,&nx,&px,&objectif,&Tol,&Maxit,bandwidth,typekernel);
  }
  if ( Maxit== (-1)) {
	Rf_error("Error: failed to find appropriate bandwidth. Try to increase dftobwitmax argument of control.par list ?");
  }
  /* result */
  UNPROTECT(1); 
  return rbandwidth;
}
