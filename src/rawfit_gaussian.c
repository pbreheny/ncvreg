#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
double crossprod(double *X, double *y, int n, int j);
double sqsum(double *x, int n, int j);
double MCP(double z, double l1, double l2, double gamma, double v);
double SCAD(double z, double l1, double l2, double gamma, double v);
double lasso(double z, double l1, double l2, double v);
double gLoss(double *r, int n);

// Memory handling, output formatting (Gaussian)
SEXP cleanupRawG(double *a, double *v, double *z, int *active, SEXP beta, SEXP loss, SEXP iter, SEXP resid) {
  Free(a);
  Free(v);
  Free(z);
  Free(active);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(res, 0, beta);
  SET_VECTOR_ELT(res, 1, loss);
  SET_VECTOR_ELT(res, 2, iter);
  SET_VECTOR_ELT(res, 3, resid);
  UNPROTECT(1);
  return(res);
}

// Coordinate descent for gaussian models
SEXP rawfit_gaussian(SEXP X_, SEXP y_, SEXP init_, SEXP r_, SEXP xtx_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_) {

  // Declarations: Outcome
  int n = length(y_);
  int p = length(X_)/n;
  SEXP res, beta, loss, iter, resid;
  PROTECT(beta = allocVector(REALSXP, p));
  double *b = REAL(beta);
  for (int j=0; j<p; j++) b[j] = 0;
  PROTECT(loss = allocVector(REALSXP, 1));
  PROTECT(iter = allocVector(INTSXP, 1));
  PROTECT(resid = allocVector(REALSXP, n));
  INTEGER(iter)[0] = 0;
  
  // Declarations
  double *X = REAL(X_);
  double *y = REAL(y_);
  double *a = Calloc(p, double); // Beta from previous iteration
  for (int j=0; j<p; j++) a[j]=REAL(init_)[j];
  const char *penalty = CHAR(STRING_ELT(penalty_, 0));
  double lam = REAL(lambda)[0];
  double eps = REAL(eps_)[0];
  int max_iter = INTEGER(max_iter_)[0];
  double gamma = REAL(gamma_)[0];
  double *m = REAL(multiplier);
  double alpha = REAL(alpha_)[0];
  int *active = Calloc(p, int);
  for (int j=0; j<p; j++) active[j] = 1*(a[j] != 0);
  double l1, l2;

  // Setup r, v, z
  double *r = REAL(resid);
  if (ISNA(REAL(r_)[0])) {
    for (int i=0; i<n; i++) r[i] = y[i];
    for (int j=0; j<p; j++) {
      for (int i=0; i<n; i++) {
        r[i] -= X[j*n+i]*a[j];
      }
    }
  } else {
    for (int i=0; i<n; i++) r[i] = REAL(r_)[i];
  }
  double *v = Calloc(p, double);
  if (ISNA(REAL(xtx_)[0])) {
    for (int j=0; j<p; j++) v[j] = sqsum(X, n, j)/n;
  } else {
    for (int j=0; j<p; j++) v[j] = REAL(xtx_)[j];
  }
  double *z = Calloc(p, double);
  double sdy = sqrt(gLoss(y, n)/n);

  // Fit
  while (INTEGER(iter)[0] < max_iter) {
    R_CheckUserInterrupt();
    while (INTEGER(iter)[0] < max_iter) {
      INTEGER(iter)[0]++;

      // Solve over the active set
      double maxChange = 0;
      for (int j=0; j<p; j++) {
        if (active[j]) {
          z[j] = crossprod(X, r, n, j)/n + v[j]*a[j];

          // Update beta_j
          l1 = lam * m[j] * alpha;
          l2 = lam * m[j] * (1-alpha);
          if (strcmp(penalty,"MCP")==0) b[j] = MCP(z[j], l1, l2, gamma, v[j]);
          if (strcmp(penalty,"SCAD")==0) b[j] = SCAD(z[j], l1, l2, gamma, v[j]);
          if (strcmp(penalty,"lasso")==0) b[j] = lasso(z[j], l1, l2, v[j]);

          // Update r
          double shift = b[j] - a[j];
          if (shift !=0) {
            for (int i=0;i<n;i++) r[i] -= shift*X[j*n+i];
            if (fabs(shift)*sqrt(v[j]) > maxChange) maxChange = fabs(shift) * sqrt(v[j]);
          }
        }
      }
        
      // Check for convergence
      for (int j=0; j<p; j++) a[j] = b[j];
      if (maxChange < eps*sdy) break;
    }

    // Scan for violations
    int violations = 0;
    for (int j=0; j<p; j++) {
      if (!active[j]) {

        z[j] = crossprod(X, r, n, j)/n;

        // Update beta_j
        l1 = lam * m[j] * alpha;
        l2 = lam * m[j] * (1-alpha);
        if (strcmp(penalty,"MCP")==0) b[j] = MCP(z[j], l1, l2, gamma, v[j]);
        if (strcmp(penalty,"SCAD")==0) b[j] = SCAD(z[j], l1, l2, gamma, v[j]);
        if (strcmp(penalty,"lasso")==0) b[j] = lasso(z[j], l1, l2, v[j]);

        // If something enters, update active set & residuals
        if (b[j] !=0) {
          active[j] = 1;
          for (int i=0; i<n; i++) r[i] -= b[j]*X[j*n+i];
          a[j] = b[j];
          violations++;
        }
      }
    }
    if (violations==0) break;
  }
  REAL(loss)[0] = gLoss(r, n);
  res = cleanupRawG(a, v, z, active, beta, loss, iter, resid);
  UNPROTECT(4);
  return(res);
}
