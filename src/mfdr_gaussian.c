#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
double crossprod(double *X, double *y, int n, int j);
SEXP getListElement(SEXP list, const char *str);
double pnorm(double x, double mu, double sigma, int lower_tail, int give_log);

SEXP mfdr_gaussian(SEXP fit) {

  // Declarations
  int n = INTEGER(getListElement(fit, "n"))[0];
  int L = ncols(getListElement(fit, "beta"));
  int p = nrows(getListElement(fit, "beta"));  
  double *b = REAL(getListElement(fit, "beta"));
  double *lambda = REAL(getListElement(fit, "lambda"));
  double *RSS = REAL(getListElement(fit, "loss"));
  double alpha = REAL(getListElement(fit, "alpha"))[0];
  double *m = REAL(getListElement(fit, "penalty.factor"));
  int S;
  double tau;
  SEXP EF;
  PROTECT(EF = allocVector(REALSXP, L));
  for (int l=0; l<L; l++) REAL(EF)[l] = 0;

  // Calculation
  for (int l=0; l<L; l++) {
    S = 1;
    for (int j=1; j<p; j++) {
      if (b[p*l+j] != 0) S += 1;
    }
    for (int j=1; j < p; j++) {
      tau = sqrt(RSS[l]/(n-S));
      REAL(EF)[l] += 2*pnorm(-sqrt(n)*lambda[l]*alpha*m[j-1]/tau, 0, 1, 1, 0);
    }
  }

  // Return
  UNPROTECT(1);
  return(EF);
}
