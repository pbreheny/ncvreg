#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
double crossprod(double *X, double *y, int n, int j);
double wcrossprod(double *X, double *y, double *w, int n, int j);
double wsqsum(double *X, double *w, int n, int j);
double sum(double *x, int n);
double MCP(double z, double l1, double l2, double gamma, double v);
double SCAD(double z, double l1, double l2, double gamma, double v);
double lasso(double z, double l1, double l2, double v);

// Memory handling, output formatting (Cox)
SEXP cleanupCox(double *a, double *r, double *h, int *e, double *eta, double *haz, double *rsk, SEXP beta, SEXP Loss, SEXP iter, SEXP Eta) {
  free(a);
  free(r);
  free(h);
  free(e);
  free(eta);
  free(haz);
  free(rsk);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(res, 0, beta);
  SET_VECTOR_ELT(res, 1, Loss);
  SET_VECTOR_ELT(res, 2, iter);
  SET_VECTOR_ELT(res, 3, Eta);
  UNPROTECT(1);
  return(res);
}

// Coordinate descent for Cox models
SEXP cdfit_cox_dh(SEXP X_, SEXP d_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_, SEXP dfmax_, SEXP user_, SEXP warn_) {

  // Lengths/dimensions
  int n = length(d_);
  int p = length(X_)/n;
  int L = length(lambda);

  // Pointers
  double *X = REAL(X_);
  double *d = REAL(d_);
  const char *penalty = CHAR(STRING_ELT(penalty_, 0));
  double *lam = REAL(lambda);
  double eps = REAL(eps_)[0];
  int max_iter = INTEGER(max_iter_)[0];
  int tot_iter = 0;
  double gamma = REAL(gamma_)[0];
  double *m = REAL(multiplier);
  double alpha = REAL(alpha_)[0];
  int dfmax = INTEGER(dfmax_)[0];
  int user = INTEGER(user_)[0];
  int warn = INTEGER(warn_)[0];

  // Outcome
  SEXP res, beta, Loss, iter, Eta;
  PROTECT(beta = allocVector(REALSXP, L*p));
  for (int j=0; j<(L*p); j++) REAL(beta)[j] = 0;
  double *b = REAL(beta);
  PROTECT(Loss = allocVector(REALSXP, L));
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;
  PROTECT(Eta = allocVector(REALSXP, L*n));
  for (int j=0; j<(L*n); j++) REAL(Eta)[j] = 0;

  // Intermediate quantities
  double *a = R_Calloc(p, double);    // Beta from previous iteration
  for (int j=0; j<p; j++) a[j] = 0;
  double *haz = R_Calloc(n, double);
  double *rsk = R_Calloc(n, double);
  double *r = R_Calloc(n, double);
  double *h = R_Calloc(n, double);
  int *e = R_Calloc(p, int);
  for (int j=0; j<p; j++) e[j] = 0;
  double *eta = R_Calloc(n, double);
  for (int i=0; i<n; i++) eta[i] = 0;
  double xwr, xwx, u, v, l1, l2, shift, si, s, nullDev;
  int lstart;

  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  rsk[n-1] = 1;
  for (int i=n-2; i>=0; i--) rsk[i] = rsk[i+1] + 1;
  nullDev = 0;
  for (int i=0; i<n; i++) nullDev -= d[i]*log(rsk[i]);
  if (user) {
    lstart = 0;
  } else {
    lstart = 1;
    REAL(Loss)[0] = nullDev;
  }

  // Path
  for (int l=lstart; l<L; l++) {
    R_CheckUserInterrupt();
    if (l != 0) {
      // Assign a
      for (int j=0; j<p; j++) a[j] = b[(l-1)*p+j];

      // Check dfmax
      int nv = 0;
      for (int j=0; j<p; j++) {
        if (a[j] != 0) nv++;
      }
      if ((nv > dfmax) | (tot_iter == max_iter)) {
        for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
        break;
      }
    }

    while (tot_iter < max_iter) {
      while (tot_iter < max_iter) {
        INTEGER(iter)[l]++;
        tot_iter++;
        REAL(Loss)[l] = 0;
        double maxChange = 0;

        // Calculate haz, risk
        for (int i=0; i<n; i++) haz[i] = exp(eta[i]);
        rsk[n-1] = haz[n-1];
        for (int i=n-2; i>=0; i--) {
          rsk[i] = rsk[i+1] + haz[i];
        }
        for (int i=0; i<n; i++) {
          REAL(Loss)[l] += d[i]*eta[i] - d[i]*log(rsk[i]);
        }

        // Approximate L
        h[0] = d[0]/rsk[0];
        for (int i=1; i<n; i++) {
          h[i] = h[i-1] + d[i]/rsk[i];
        }
        for (int i=0; i<n; i++) {
          h[i] = h[i]*haz[i];
          s = d[i] - h[i];
          if (h[i]==0) r[i]=0;
          else r[i] = s/h[i];
        }

        // Check for saturation
        if (REAL(Loss)[l]/nullDev < .01) {
          if (warn) warning("Model saturated; exiting...");
          for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
          tot_iter = max_iter;
          break;
        }

        // Covariates
        for (int j=0; j<p; j++) {
          if (e[j]) {

            // Calculate u, v
            xwr = wcrossprod(X, r, h, n, j);
            xwx = wsqsum(X, h, n, j);
            u = xwr/n + (xwx/n)*a[j];
            v = xwx/n;

            // Update b_j
            l1 = lam[l] * m[j] * alpha;
            l2 = lam[l] * m[j] * (1-alpha);
            if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(u, l1, l2, gamma, v);
            if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(u, l1, l2, gamma, v);
            if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(u, l1, l2, v);

            // Update r
            shift = b[l*p+j] - a[j];
            if (shift !=0) {
              for (int i=0;i<n;i++) {
                si = shift*X[j*n+i];
                r[i] -= si;
                eta[i] += si;
              }
              if (fabs(shift)*sqrt(v) > maxChange) maxChange = fabs(shift)*sqrt(v);
            }
          }
        }

        // Check for convergence
        for (int j=0; j<p; j++) a[j] = b[l*p+j];
        if (maxChange < eps) break;
      }

      // Scan for violations
      int violations = 0;
      for (int j=0; j<p; j++) {
        if (e[j]==0) {
          xwr = wcrossprod(X, r, h, n, j)/n;
          l1 = lam[l] * m[j] * alpha;
          if (fabs(xwr) > l1) {
            e[j] = 1;
            violations++;
          }
        }
      }
      if (violations==0) {
        for (int i=0; i<n; i++) {
          REAL(Eta)[l*n+i] = eta[i];
        }
        break;
      }
    }
  }
  res = cleanupCox(a, r, h, e, eta, haz, rsk, beta, Loss, iter, Eta);
  UNPROTECT(4);
  return(res);
}
