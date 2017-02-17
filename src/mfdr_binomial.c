#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
double crossprod(double *X, double *y, int n, int j);
double wsqsum(double *X, double *w, int n, int j);
SEXP getListElement(SEXP list, const char *str);
double p_binomial(double eta);
double pnorm(double x, double mu, double sigma, int lower_tail, int give_log);

SEXP mfdr_binomial(SEXP fit) {

  // Declarations
  int L = length(getListElement(fit, "lambda"));
  int n = INTEGER(getListElement(fit, "n"))[0];
  int p = ncols(getListElement(fit, "X"));
  double *X = REAL(getListElement(fit, "X"));
  double *Eta = REAL(getListElement(fit, "Eta"));
  double *lambda = REAL(getListElement(fit, "lambda"));
  double alpha = REAL(getListElement(fit, "alpha"))[0];
  double *m = REAL(getListElement(fit, "penalty.factor"));
  double pi, tau;
  double *w = Calloc(n, double);
  SEXP EF;
  PROTECT(EF = allocVector(REALSXP, L));
  for (int l=0; l<L; l++) REAL(EF)[l] = 0;

  // Calculation
  for (int l=0; l<L; l++) {
    for (int i=0; i<n; i++) {
      pi = p_binomial(Eta[n*l+i]);
      w[i] = pi*(1-pi);
    }
    for (int j=0; j<p; j++) {
      tau = sqrt(wsqsum(X, w, n, j)/n);
      REAL(EF)[l] += 2*pnorm(-sqrt(n)*lambda[l]*alpha*m[j]/tau, 0, 1, 1, 0);
    }
  }
  // Rprintf("L: %d\n", L);

  // Return
  Free(w);
  UNPROTECT(1);
  return(EF);
}
