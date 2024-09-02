#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
double crossprod(double *X, double *y, int n, int j);
SEXP getListElement(SEXP list, const char *str);
double wsqsum(double *X, double *w, int n, int j);
double pnorm(double x, double mu, double sigma, int lower_tail, int give_log);

SEXP mfdr_cox(SEXP fit) {

  // Declarations
  int L = length(getListElement(fit, "lambda"));
  int n = INTEGER(getListElement(fit, "n"))[0];
  int p = ncols(getListElement(fit, "X"));
  double *X = REAL(getListElement(fit, "X"));
  double *d = REAL(getListElement(fit, "fail"));
  double *Eta = REAL(getListElement(fit, "linear.predictors"));
  double *lambda = REAL(getListElement(fit, "lambda"));
  double alpha = REAL(getListElement(fit, "alpha"))[0];
  double *m = REAL(getListElement(fit, "penalty.factor"));
  double tau;
  double *w = R_Calloc(n, double);
  double *haz = R_Calloc(n, double);
  double *rsk = R_Calloc(n, double);
  SEXP EF;
  PROTECT(EF = allocVector(REALSXP, L));
  for (int l=0; l<L; l++) REAL(EF)[l] = 0;

  // Calculation
  for (int l=0; l<L; l++) {
    for (int i=0; i<n; i++) haz[i] = exp(Eta[n*l+i]);
    rsk[n-1] = haz[n-1];
    for (int i=n-2; i>=0; i--) rsk[i] = rsk[i+1] + haz[i];
    for (int j=0; j<n; j++) {
      w[j] = 0;
      for (int i=0; i <= j; i++) {
        w[j] += d[i]*haz[j]/rsk[i]*(1-haz[j]/rsk[i]);
      }
    }
    for (int j=0; j<p; j++) {
      tau = sqrt(wsqsum(X, w, n, j)/n);
      REAL(EF)[l] += 2*pnorm(-sqrt(n)*lambda[l]*alpha*m[j]/tau, 0, 1, 1, 0);
    }
  }

  // Return
  free(w);
  free(haz);
  free(rsk);
  UNPROTECT(1);
  return(EF);
}
