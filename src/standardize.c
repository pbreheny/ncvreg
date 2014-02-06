#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>

SEXP standardize(SEXP X_) {
  // Declarations
  int n = nrows(X_);
  int p = ncols(X_);
  SEXP XX_, c_, s_;
  PROTECT(XX_ = allocMatrix(REALSXP, n, p));
  PROTECT(c_ = allocVector(REALSXP, p));
  PROTECT(s_ = allocVector(REALSXP, p));
  double *X = REAL(X_);
  double *XX = REAL(XX_);
  double *c = REAL(c_);
  double *s = REAL(s_);

  for (int j=0; j<p; j++) {
    // Center
    c[j] = 0;
    for (int i=0; i<n; i++) {
      c[j] += X[j*n+i];
    }
    c[j] = c[j] / n;
    for (int i=0; i<n; i++) XX[j*n+i] = X[j*n+i] - c[j];

    // Scale
    s[j] = 0;
    for (int i=0; i<n; i++) {
      s[j] += pow(XX[j*n+i], 2);
    }
    s[j] = sqrt(s[j]/n);
    for (int i=0; i<n; i++) XX[j*n+i] = XX[j*n+i]/s[j];
  }

  // Return list
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(res, 0, XX_);
  SET_VECTOR_ELT(res, 1, c_);
  SET_VECTOR_ELT(res, 2, s_);
  UNPROTECT(4);
  return(res);
}
