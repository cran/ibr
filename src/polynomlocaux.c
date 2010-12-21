#include <math.h>
/* Regression noyau gaussien */
void regpolg(double *x, int *nx,double *y, double *bw, double *valx, int *nvalx, double *regx,  double *df, double *deriv)
{
  int i, j ; 

  double S0, S1, S2, w, T0, T1, wii;
  /* initialisation */      
  w = 0.0;
  *df=0.0;
  for(i = 0; i < *nvalx; i++)
    regx[i] = 0.0;
  for(i = 0; i < *nvalx; i++) {
    wii=0.0;
    S0 = 0.0;
    S1 = 0.0;
    S2 = 0.0;
    T1 = 0.0;
    T0 = 0.0;
    /* pour la i eme valeur de la grille valx :*/
    /* boucle sur les valeurs observees (indice j)*/
    for(j = 0; j < *nx; j++) {
      /* poids */
      w= exp(-0.5*(pow((valx[i]-x[j])/(*bw) ,2))) /sqrt(2*3.14159265358979);
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
      *df=*df+wii/S0;
    }
  }
}

void regpole(double *x, int *nx,double *y, double *bw, double *valx, int *nvalx, double *regx,  double *df, double *deriv)
{
  int i, j ; 

  double S0, S1, S2, w, T0, T1, xc, wii;
  /* initialisation */      
  w = 0.0;
  xc = 0.0;
  *df=0.0;
  for(i = 0; i < *nvalx; i++)
    regx[i] = 0.0;
  for(i = 0; i < *nvalx; i++) {
    wii=0.0;
    S0 = 0.0;
    S1 = 0.0;
    S2 = 0.0;
    T1 = 0.0;
    T0 = 0.0;
    /* pour la i eme valeur de la grille valx :*/
    /* boucle sur les valeurs observees (indice j)*/
    for(j = 0; j < *nx; j++) {
      /* poids */
      xc = pow((valx[i]-x[j])/(*bw) ,2);
      if (xc<= 1.0) {
	w=3.0 / 4.0 * (1.0 - xc) ;
	if (i==j) wii=w;
	/* sommes */ 
	S0 = S0 + w;
	S1 = S1 + w * (x[j] - valx[i]);
	S2 = S2 + w * pow((x[j] - valx[i]),2);
	/* regression */ 
	T0 = T0+w*y[j];
	T1 = T1+ (x[j] - valx[i]) * w * y[j];
      }
    }
    if (S0>0) {
      regx[i]= (S2 * T0 - S1 * T1)/(S0 * S2 - pow(S1,2));
      deriv[i]=(- S1 * T0 + S0 * T1)/(S0 * S2 - pow(S1,2));
      *df=*df+wii/S0;
    }
  }
}



void regpolq(double *x, int *nx,double *y, double *bw, double *valx, int *nvalx, double *regx,  double *df, double *deriv)
{
  int i, j ; 

  double S0, S1, S2, w, T0, T1, xc, wii;
  /* initialisation */      
  w = 0.0;
  xc = 0.0;
  *df=0.0;
  for(i = 0; i < *nvalx; i++)
    regx[i] = 0.0;
  for(i = 0; i < *nvalx; i++) {
    wii=0.0;
    S0 = 0.0;
    S1 = 0.0;
    S2 = 0.0;
    T1 = 0.0;
    T0 = 0.0;
    /* pour la i eme valeur de la grille valx :*/
    /* boucle sur les valeurs observees (indice j)*/
    for(j = 0; j < *nx; j++) {
      /* poids */
      xc = pow((valx[i]-x[j])/(*bw) ,2);
      if (xc<= 1.0) {
	w=15.0 / 16.0 * pow(pow((1.0 - xc),2),2) ;
	if (i==j) wii=w;
	/* sommes */ 
	S0 = S0 + w;
	S1 = S1 + w * (x[j] - valx[i]);
	S2 = S2 + w * pow((x[j] - valx[i]),2);
	/* regression */ 
	T0 = T0+w*y[j];
	T1 = T1+ (x[j] - valx[i]) * w * y[j];
      }
    }
    if (S0>0) {
      regx[i] = (S2 * T0 - S1 * T1)/(S0 * S2 - pow(S1,2));
      deriv[i] = (- S1 * T0 + S0 * T1)/(S0*S2 - pow(S1,2));
      *df=*df+wii/S0;
    }
  }
}



void regpolu(double *x, int *nx,double *y, double *bw, double *valx, int *nvalx, double *regx,  double *df, double *deriv)
{
  int i, j ; 

  double S0, S1, S2, w, T0, T1, xc, wii;
  /* initialisation */      
  w = 0.0;
  xc = 0.0;
  *df=0.0;
  for(i = 0; i < *nvalx; i++)
    regx[i] = 0.0;
  for(i = 0; i < *nvalx; i++) {
    wii=0.0;
    S0 = 0.0;
    S1 = 0.0;
    S2 = 0.0;
    T1 = 0.0;
    T0 = 0.0;
    /* pour la i eme valeur de la grille valx :*/
    /* boucle sur les valeurs observees (indice j)*/
    for(j = 0; j < *nx; j++) {
      /* poids */
      xc = fabs((valx[i]-x[j])/(*bw));
      if (xc<= 1.0) {
	w=0.5 ;
	if (i==j) wii=w;
	/* sommes */ 
	S0 = S0 + w;
	S1 = S1 + w * (x[j] - valx[i]);
	S2 = S2 + w * pow((x[j] - valx[i]),2);
	/* regression */ 
	T0 = T0+w*y[j];
	T1 = T1+ (x[j] - valx[i]) * w * y[j];
      }
    }
    if (S0>0) {
      regx[i] = (S2 * T0 - S1 * T1)/(S0 * S2 - pow(S1,2));
      deriv[i] = (- S1 * T0 + S0 * T1)/(S0*S2 - pow(S1,2));
      *df=*df+wii/S0;
    }
  }
}

/***********************************************************/
/***********************************************************/
/***********************************************************/
/* Regression noyau gaussien */
void regpolgcv(double *x, int *nx,double *y, double *bw, int *nbw,  double *effold, int *neffold, double *sse, double *sap )
{
  int i, j , k, h ; 
  
  double S0, S1, S2, w, T0, T1,  regx;
  /* initialisation */      
  w = 0.0;
  /* boucle sur les fenetres*/
  for (h=0; h < *nbw; h++) {
    sse[h]=0.0;
    sap[h]=0.0;
    /* boucle sur les fold */
    for (k=0; k < *neffold; k++) {
      /* boucle sur la partie test */
      for(i = effold[k]; i < effold[k+1]; i++) {
	S0 = 0.0;
	S1 = 0.0;
	S2 = 0.0;
	T1 = 0.0;
	T0 = 0.0;
	/* pour la i eme valeur de la grille x :*/
	/* boucle sur les valeurs observees (indice j)*/
	for(j = 0; j < *nx; j++) {
	  /* si pas en test */
	  if ((j>=effold[k+1])||(j<effold[k])) {
	    /* poids */
	    w= exp(-0.5*(pow((x[i]-x[j])/(bw[h]) ,2))) /sqrt(2*3.14159265358979);
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
}

/***********************************************************/
void regpolecv(double *x, int *nx,double *y, double *bw, int *nbw, double *effold, int *neffold, double *sse, double *sap )
{
  int i, j , k ,h; 

  double S0, S1, S2, w, T0, T1, xc, regx;
  for (h=0; h < *nbw; h++) {
    sse[h]=0.0;
    sap[h]=0.0;
    /* boucle sur les fold */
    for (k=0; k < *neffold; k++) {
      /* boucle sur la partie test */
      for(i = effold[k]; i < effold[k+1]; i++) {
	w = 0.0;
	S0 = 0.0;
	S1 = 0.0;
	S2 = 0.0;
	T1 = 0.0;
	T0 = 0.0;
	/* pour la i eme valeur de la grille x :*/
	/* boucle sur les valeurs observees (indice j)*/
	for(j = 0; j < *nx; j++) {
	  /* si pas en test */
	  if ((j>=effold[k+1])||(j<effold[k])) {
	    /* poids */
	    xc = pow((x[i]-x[j])/(bw[h]) ,2);
	    if (xc<= 1.0) {
	      w=3.0 / 4.0 * (1.0 - xc) ;
	      /* sommes */ 
	      S0 = S0 + w;
	      S1 = S1 + w * (x[j] - x[i]);
	      S2 = S2 + w * pow((x[j] - x[i]),2);
	      /* regression */ 
	      T0 = T0+w*y[j];
	      T1 = T1+ (x[j] - x[i]) * w * y[j];
	    }
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
}

/***********************************************************/
void regpolqcv(double *x, int *nx,double *y, double *bw, int *nbw, double *effold, int *neffold, double *sse, double *sap )
{
  int i, j , k ,h; 

  double S0, S1, S2, w, T0, T1, xc, regx;
  for (h=0; h < *nbw; h++) {
    sse[h]=0.0;
    sap[h]=0.0;
    /* boucle sur les fold */
    for (k=0; k < *neffold; k++) {
      /* boucle sur la partie test */
      for(i = effold[k]; i < effold[k+1]; i++) {
	w = 0.0;
	S0 = 0.0;
	S1 = 0.0;
	S2 = 0.0;
	T1 = 0.0;
	T0 = 0.0;
	/* pour la i eme valeur de la grille x :*/
	/* boucle sur les valeurs observees (indice j)*/
	for(j = 0; j < *nx; j++) {
	  /* si pas en test */
	  if ((j>=effold[k+1])||(j<effold[k])) {
	    /* poids */
	    xc = pow((x[i]-x[j])/(bw[h]) ,2);
	    if (xc<= 1.0) {
	      w=15.0 / 16.0 * pow(pow((1.0 - xc),2),2) ;
	      /* sommes */ 
	      S0 = S0 + w;
	      S1 = S1 + w * (x[j] - x[i]);
	      S2 = S2 + w * pow((x[j] - x[i]),2);
	      /* regression */ 
	      T0 = T0+w*y[j];
	      T1 = T1+ (x[j] - x[i]) * w * y[j];
	    }
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
}

/***********************************************************/
void regpolucv(double *x, int *nx,double *y, double *bw, int *nbw, double *effold, int *neffold, double *sse, double *sap )
{
  int i, j , k ,h; 

  double S0, S1, S2, w, T0, T1, xc, regx;
  for (h=0; h < *nbw; h++) {
    sse[h]=0.0;
    sap[h]=0.0;
    /* boucle sur les fold */
    for (k=0; k < *neffold; k++) {
      /* boucle sur la partie test */
      for(i = effold[k]; i < effold[k+1]; i++) {
	w = 0.0;
	S0 = 0.0;
	S1 = 0.0;
	S2 = 0.0;
	T1 = 0.0;
	T0 = 0.0;
	/* pour la i eme valeur de la grille x :*/
	/* boucle sur les valeurs observees (indice j)*/
	for(j = 0; j < *nx; j++) {
	  /* si pas en test */
	  if ((j>=effold[k+1])||(j<effold[k])) {
	    /* poids */
	    xc = pow((x[i]-x[j])/(bw[h]) ,2);
	    if (xc<= 1.0) {
	      w=0.5 ;
	      /* sommes */ 
	      S0 = S0 + w;
	      S1 = S1 + w * (x[j] - x[i]);
	      S2 = S2 + w * pow((x[j] - x[i]),2);
	      /* regression */ 
	      T0 = T0+w*y[j];
	      T1 = T1+ (x[j] - x[i]) * w * y[j];
	    }
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
}
