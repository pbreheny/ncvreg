#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
double crossprod(double *X, double *y, int n, int j);

SEXP maxprod(SEXP X_, SEXP y_, SEXP v_, SEXP m_) {

  // Declarations
  int n = nrows(X_);
  int p = length(v_);
  SEXP z;
  PROTECT(z = allocVector(REALSXP, 1));
  REAL(z)[0] = 0;
  double zz;
  double *X = REAL(X_);
  double *y = REAL(y_);
  double *m = REAL(m_);
  int *v = INTEGER(v_);

  for (int j=0; j<p; j++) {
    zz = crossprod(X, y, n, v[j]-1) / m[v[j]-1];
    if (fabs(zz) > REAL(z)[0]) REAL(z)[0] = fabs(zz);
  }

  // Return list
  UNPROTECT(1);
  return(z);
}
