#include "kernel.h"
#include <math.h>
#include <R.h>
#include <Rinternals.h>
const double MY_SQRT2PI=2.5066282746310005024157652848111;
static R_INLINE double kernelgauss(double x, double y, double bw)
{
  double res;
  res = exp(-0.5*(pow((x-y)/(bw) ,2))) / MY_SQRT2PI;
  return res;
}
static R_INLINE double kernelepane(double x, double y, double bw)
{
  double xc, res;
  xc = pow((x-y)/bw ,2);
  if (xc<= 1.0) {
    res = 3.0 / 4.0 * (1.0 - xc);
    return res;
  }
  else {
    return 0.0;
  }
}
static R_INLINE double kernelquad(double x, double y, double bw)
{
  double xc, res;
  xc = pow((x-y)/bw ,2);
  if (xc<= 1.0) {
    res = 15.0 / 16.0 * pow(1.0 - xc, 2);
    return res;
  }
  else {
    return 0.0;
  }
}
static R_INLINE double kernelunif(double x, double y, double bw)
{
  double xc, res;
  xc = fabs((x-y)/bw);
  if (xc<= 1.0) {
    res = 0.5; 
    return res;
  }
  else {
    return 0.0;
  }
}

double poidskernel(double x, double y, double bw, int typekernel){
  double w;
  switch(typekernel) {
  case 1:
    /* gaussian */
    w= kernelgauss(x, y, bw);
    break;
  case 2:
    /* epanechnikov */
    w= kernelepane(x, y, bw);
    break;
  case 3:
    /* quadratic */
    w= kernelquad(x, y, bw);
    break;
  case 4:
    /* uniform */
    w= kernelunif(x, y, bw);
    break;
  default:
    w=0.0;
    break;
  }
  return w;
}
