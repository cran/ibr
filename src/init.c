#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>  // optional
#include "ibr.h"
#include "product.h"
#include "regression.h"
#include "dssplines.h"
#include "smoother.h"
#include "kernel.h"
#include "polynomlocaux.h"
#include "regression.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef R_CallDef[] = {
  CALLDEF(Amatrix, 1),
  CALLDEF(Smatrix, 1),
  CALLDEF(Kmatrix, 4),
  CALLDEF(choosebw, 5),
  CALLDEF(product, 3),
  CALLDEF(productandcrit, 5),
  CALLDEF(polynom, 2),
  CALLDEF(semikerlog, 4),
  CALLDEF(semikerpow, 4),
  CALLDEF(npreg, 5),
  CALLDEF(npregcv, 5),
  CALLDEF(npregpol, 5),
  CALLDEF(npregpolcv, 5),
   {NULL, NULL, 0}
};

void
attribute_visible  // optional
R_init_ibr(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

