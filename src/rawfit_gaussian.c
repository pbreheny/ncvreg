#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
double crossprod(double *X, double *y, int n, int j);
double sum(double *x, int n);
double MCP(double z, double l1, double l2, double gamma, double v);
double SCAD(double z, double l1, double l2, double gamma, double v);
double lasso(double z, double l1, double l2, double v);

// Memory handling, output formatting (Gaussian)
SEXP cleanupG(double *a, double *r, int *e1, int *e2, double *z, SEXP beta, SEXP loss, SEXP iter) {
  Free(a);
  Free(r);
  Free(e1);
  Free(e2);
  Free(z);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(res, 0, beta);
  SET_VECTOR_ELT(res, 1, loss);
  SET_VECTOR_ELT(res, 2, iter);
  UNPROTECT(1);
  return(res);
}

// Gaussian loss
double gLoss(double *r, int n) {
  double l = 0;
  for (int i=0;i<n;i++) l = l + pow(r[i],2);
  return(l);
}

// Coordinate descent for gaussian models
SEXP rawfit_gaussian(SEXP X_, SEXP y_, SEXP init_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_) {

  // Declarations: Outcome
  int n = length(y_);
  int p = length(X_)/n;
  SEXP res, beta, loss, iter;
  PROTECT(beta = allocVector(REALSXP, p));
  double *b = REAL(beta);
  for (int j=0; j<(L*p); j++) b[j] = 0;
  PROTECT(loss = allocVector(REALSXP, 1));
  PROTECT(iter = allocVector(INTSXP, 1));
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
  double l1, l2;

  // Setup r, v, z  
  double *r = Calloc(n, double);
  for (int i=0; i<n; i++) r[i] = y[i];
  for (int j=0; j<p; j++) {
    for (int i=0; i<n; i++) {
      r[i] -= X[j*n+i]*a[j];
    }
  }
  double *v = Calloc(p, double);
  for (int j=0; j<p; j++) v[j] = sqsum(X, n, j);
  double *z = Calloc(p, double);
  double sdy = sqrt(rss/n);

  // Fit
  while (INTEGER(iter)[0] < max_iter) {
    R_CheckUserInterrupt();
    while (INTEGER(iter)[0] < max_iter) {
      INTEGER(iter)[0]++;      

      // Solve over the active set
      double maxChange = 0;
      for (int j=0; j<p; j++) {
        if (active[j]) {
          z[j] = crossprod(X, r, n, j)/v[j] + a[j];

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
            if (fabs(shift) > maxChange) maxChange = fabs(shift);
          }

        }
        
        // Check for convergence
        for (int j=0; j<p; j++) a[j] = b[l*p+j];
        if (maxChange < eps*sdy) break;

        // Scan for violations in strong set
        int violations = 0;
        for (int j=0; j<p; j++) {
          if (e1[j]==0 && e2[j]==1) {

            z[j] = crossprod(X, r, n, j)/n;

            // Update beta_j
            l1 = lam[l] * m[j] * alpha;
            l2 = lam[l] * m[j] * (1-alpha);
            if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(z[j], l1, l2, gamma, 1);
            if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(z[j], l1, l2, gamma, 1);
            if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(z[j], l1, l2, 1);

            // If something enters the eligible set, update eligible set & residuals
            if (b[l*p+j] !=0) {
              e1[j] = e2[j] = 1;
              for (int i=0; i<n; i++) r[i] -= b[l*p+j]*X[j*n+i];
              a[j] = b[l*p+j];
              violations++;
            }
          }
        }
        if (violations==0) break;
      }

      // Scan for violations in rest
      int violations = 0;
      for (int j=0; j<p; j++) {
        if (e2[j]==0) {

          z[j] = crossprod(X, r, n, j)/n;

          // Update beta_j
          l1 = lam[l] * m[j] * alpha;
          l2 = lam[l] * m[j] * (1-alpha);
          if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(z[j], l1, l2, gamma, 1);
          if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(z[j], l1, l2, gamma, 1);
          if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(z[j], l1, l2, 1);

          // If something enters the eligible set, update eligible set & residuals
          if (b[l*p+j] !=0) {
            e1[j] = e2[j] = 1;
            for (int i=0; i<n; i++) r[i] -= b[l*p+j]*X[j*n+i];
            a[j] = b[l*p+j];
            violations++;
          }
        }
      }

      if (violations==0) {
        break;
      }
    }
    REAL(loss)[l] = gLoss(r, n);
  }
  res = cleanupG(a, r, e1, e2, z, beta, loss, iter);
  UNPROTECT(3);
  return(res);
}
