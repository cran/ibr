#include <math.h>
/* Regression noyau gaussien */
void regg(double *x, int *nx,double *y, double *bw, double *valx, int *nvalx, double *regx, double *df)
    {
      int i, j ; 

      double some, w, wii;
      /* initialisation */      
      w = 0.0;
      *df=0.0;
      for(i = 0; i < *nvalx; i++)
        regx[i] = 0.0;
      for(i = 0; i < *nvalx; i++) {
	wii=0.0;
	some = 0.0;
      /* pour la i eme valeur de la grille valx :*/
      /* boucle sur les valeurs observees (indice j)*/
	for(j = 0; j < *nx; j++) {
      /* poids */
	  w= exp(-0.5*(pow((valx[i]-x[j])/(*bw) ,2))) /sqrt(2*3.14159265358979);
	  if (i==j) wii=w;
	  some=some+w;
      /* regression */ 
	  regx[i]=regx[i]+w*y[j];
	}
	if (some>0) {
	  regx[i]=regx[i]/some;
	  *df=*df+wii/some;
	}
      }
    }

void rege(double *x, int *nx,double *y, double *bw, double *valx, int *nvalx, double *regx, double *df)
    {
      int i, j ; 

      double some, w, xc, wii;
      /* initialisation */      
      w = 0.0;
      xc = 0.0;
      *df=0.0;
      for(i = 0; i < *nvalx; i++)
        regx[i] = 0.0;
      for(i = 0; i < *nvalx; i++) {
	wii=0.0;
	some = 0.0;
      /* pour la i eme valeur de la grille valx :*/
      /* boucle sur les valeurs observees (indice j)*/
	for(j = 0; j < *nx; j++) {
      /* poids */
	  xc = pow((valx[i]-x[j])/(*bw) ,2);
	  if (xc<= 1.0) {
	    w=3.0 / 4.0 * (1.0 - xc) ;
	    if (i==j) wii=w;
	    some=some+w; 
	    /* regression */ 
	    regx[i]=regx[i]+w*y[j];
	  }
	}
	if (some>0) {
	  regx[i]=regx[i]/some;
	  *df=*df+wii/some;
	}
      }
    }



void regq(double *x, int *nx,double *y, double *bw, double *valx, int *nvalx, double *regx, double *df)
    {
      int i, j ; 

      double some, w, xc, wii;
      /* initialisation */      
      w = 0.0;
      xc = 0.0;
      *df=0.0;
      for(i = 0; i < *nvalx; i++)
        regx[i] = 0.0;
      for(i = 0; i < *nvalx; i++) {
	wii=0.0;
	some = 0.0;
      /* pour la i eme valeur de la grille valx :*/
      /* boucle sur les valeurs observees (indice j)*/
	for(j = 0; j < *nx; j++) {
      /* poids */
	  xc = pow((valx[i]-x[j])/(*bw) ,2);
	  if (xc<= 1.0) {
	    w=15.0 / 16.0 * pow(pow((1.0 - xc),2),2) ;
	    if (i==j) wii=w;
	    /* regression */ 
	    some=some+w; 
	    regx[i]=regx[i]+w*y[j];
	  }
	}
	if (some>0) {
	  regx[i]=regx[i]/some;
	  *df=*df+wii/some;
	}
      }
    }

void regu(double *x, int *nx,double *y, double *bw, double *valx, int *nvalx, double *regx, double *df)
    {
      int i, j ; 

      double some, w, xc, wii;
      /* initialisation */      
      w = 0.0;
      xc = 0.0;
      *df=0.0;
      for(i = 0; i < *nvalx; i++)
        regx[i] = 0.0;
      for(i = 0; i < *nvalx; i++) {
	wii=0.0;
	some = 0.0;
      /* pour la i eme valeur de la grille valx :*/
      /* boucle sur les valeurs observees (indice j)*/
	for(j = 0; j < *nx; j++) {
      /* poids */
	  xc = fabs((valx[i]-x[j])/(*bw));
	  if (xc<= 1.0) {
	    w=0.5 ;
	    if (i==j) wii=w;
	    /* regression */ 
	    some=some+w; 
	    regx[i]=regx[i]+w*y[j];
	  }
	}
	if (some>0) {
	  regx[i]=regx[i]/some;
	  *df=*df+wii/some;
	}
      }
    }

/***********************************************************/
/***********************************************************/
/***********************************************************/
void reggcv(double *x, int *nx,double *y, double *bw, int *nbw, double *effold, int *neffold, double *sse, double *sap )
{
  int i, j, k, h ; 

  double some, w, regx;
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
	regx = 0.0;
	some = 0.0;
	/* pour la i eme valeur de la valeur en test :*/
	/* boucle sur les valeurs en apprentissage (indice j)*/
	for(j = 0; j < *nx; j++ ) {
	  /* si pas en test */
	  if ((j>=effold[k+1])||(j<effold[k])) {
	    /* alors calculs */
	    /* poids */
	    w= exp(-0.5*(pow((x[i]-x[j])/(bw[h]) ,2))) /sqrt(2*3.14159265358979);
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
}

/***********************************************************/
void regecv(double *x, int *nx,double *y, double *bw, int *nbw, double *effold, int *neffold, double *sse, double *sap )
{
  int i, j, k, h ; 

  double some, w, xc, regx;
  /* boucle sur les fenetres*/
  for (h=0; h < *nbw; h++) {
    sse[h]=0.0;
    sap[h]=0.0;
    /* boucle sur les fold */
    for (k=0; k < *neffold; k++) {
      /* boucle sur la partie test */
      for(i = effold[k]; i < effold[k+1]; i++) {
	/* initialisation */      
	w = 0.0;
	regx = 0.0;
	some = 0.0;
	/* pour la i eme valeur de la valeur en test :*/
	/* boucle sur les valeurs en apprentissage (indice j)*/
	for(j = 0; j < *nx; j++ ) {
	  /* si pas en test */
	  if ((j>=effold[k+1])||(j<effold[k])) {
	    /* alors calculs */
	    /* poids */
	    xc = pow((x[i]-x[j])/(bw[h]) ,2);
	    if (xc<= 1.0) {
	    /* poids */
	    w= 3.0 / 4.0 * (1.0 - xc) ;
	    some=some+w;
	    /* regression */ 
	    regx=regx+w*y[j];
	    }
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
}
/***********************************************************/
void regqcv(double *x, int *nx,double *y, double *bw, int *nbw, double *effold, int *neffold, double *sse, double *sap )
{
  int i, j, k, h ; 

  double some, w, xc, regx;
  /* initialisation */      
  w = 0.0;
  /* boucle sur les fenetres*/
  for (h=0; h < *nbw; h++) {
    /* boucle sur les fold */
    for (k=0; k < *neffold; k++) {
      /* boucle sur la partie test */
      for(i = effold[k]; i < effold[k+1]; i++) {
	regx = 0.0;
	some = 0.0;
	/* pour la i eme valeur de la valeur en test :*/
	/* boucle sur les valeurs en apprentissage (indice j)*/
	for(j = 0; j < *nx; j++ ) {
	  /* si pas en test */
	  if ((j>=effold[k+1])||(j<effold[k])) {
	    /* alors calculs */
	    /* poids */
	    xc = pow((x[i]-x[j])/(bw[h]) ,2);
	    if (xc<= 1.0) {
	    /* poids */
	    w= 15.0 / 16.0 * pow(pow((1.0 - xc),2),2) ;
	    some=some+w;
	    /* regression */ 
	    regx=regx+w*y[j];
	    }
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
}
/***********************************************************/
void regucv(double *x, int *nx,double *y, double *bw, int *nbw, double *effold, int *neffold, double *sse, double *sap )
{
  int i, j, k, h ; 

  double some, w, xc, regx;
  /* initialisation */      
  w = 0.0;
  /* boucle sur les fenetres*/
  for (h=0; h < *nbw; h++) {
    /* boucle sur les fold */
    for (k=0; k < *neffold; k++) {
      /* boucle sur la partie test */
      for(i = effold[k]; i < effold[k+1]; i++) {
	regx = 0.0;
	some = 0.0;
	/* pour la i eme valeur de la valeur en test :*/
	/* boucle sur les valeurs en apprentissage (indice j)*/
	for(j = 0; j < *nx; j++ ) {
	  /* si pas en test */
	  if ((j>=effold[k+1])||(j<effold[k])) {
	    /* alors calculs */
	    /* poids */
	    xc = pow((x[i]-x[j])/(bw[h]) ,2);
	    if (xc<= 1.0) {
	    /* poids */
	    w= 0.5 ;
	    some=some+w;
	    /* regression */ 
	    regx=regx+w*y[j];
	    }
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
}
