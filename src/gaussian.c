#include <math.h>
#include <string.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R.h>
#include <R_ext/Applic.h>
double crossprod(double *X, double *y, int n, int j);
double sum(double *x, int n);
double MCP(double z, double l1, double l2, double gamma, double v);
double SCAD(double z, double l1, double l2, double gamma, double v);
double lasso(double z, double l1, double l2, double v);
double g_loss(double *r, int n);

// Memory handling, output formatting (Gaussian)
SEXP cleanup_gaussian(double *a, double *r, int *e1, int *e2, double *z, SEXP beta, SEXP loss, SEXP Eta, SEXP iter) {
  Free(a);
  Free(r);
  Free(e1);
  Free(e2);
  Free(z);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(res, 0, beta);
  SET_VECTOR_ELT(res, 1, loss);
  SET_VECTOR_ELT(res, 2, Eta);
  SET_VECTOR_ELT(res, 3, iter);
  UNPROTECT(1);
  return(res);
}

// Coordinate descent for gaussian models
SEXP cdfit_gaussian(SEXP X_, SEXP y_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_, SEXP dfmax_, SEXP user_) {

  // Lengths/dimensions
  int n = length(y_);
  int p = length(X_)/n;
  int L = length(lambda);

  // Pointers
  double *X = REAL(X_);
  double *y = REAL(y_);
  const char *penalty = CHAR(STRING_ELT(penalty_, 0));
  double *lam = REAL(lambda);
  double eps = REAL(eps_)[0];  
  double gamma = REAL(gamma_)[0];
  double *m = REAL(multiplier);
  double alpha = REAL(alpha_)[0];
  int max_iter = INTEGER(max_iter_)[0];
  int tot_iter = 0;
  int dfmax = INTEGER(dfmax_)[0];
  int user = INTEGER(user_)[0];

  // Outcome
  SEXP res, beta, loss, Eta, iter;
  PROTECT(beta = allocVector(REALSXP, L*p));
  for (int j=0; j<(L*p); j++) REAL(beta)[j] = 0;
  double *b = REAL(beta);
  PROTECT(loss = allocVector(REALSXP, L));
  for (int i=0; i<L; i++) REAL(loss)[i] = 0;
  PROTECT(Eta = allocVector(REALSXP, L*n));
  for (int j=0; j<(L*n); j++) REAL(Eta)[j] = 0;
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;

  // Intermediate quantities
  double *a = Calloc(p, double);
  for (int j=0; j<p; j++) a[j]=0;
  double *r = Calloc(n, double);
  for (int i=0; i<n; i++) r[i] = y[i];
  double *z = Calloc(p, double);
  for (int j=0; j<p; j++) z[j] = crossprod(X, r, n, j)/n;
  int *e1 = Calloc(p, int);
  for (int j=0; j<p; j++) e1[j] = 0;
  int *e2 = Calloc(p, int);
  for (int j=0; j<p; j++) e2[j] = 0;
  double cutoff, l1, l2;
  int lstart;

  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  double rss = g_loss(r,n);
  if (user) {
    lstart = 0;
  } else {
    REAL(loss)[0] = rss;
    lstart = 1;
  }
  double sdy = sqrt(rss/n);

  // Path
  for (int l=lstart;l<L;l++) {
    R_CheckUserInterrupt();
    if (l != 0) {
      // Assign a
      for (int j=0;j<p;j++) a[j] = b[(l-1)*p+j];

      // Check dfmax
      int nv = 0;
      for (int j=0; j<p; j++) {
        if (a[j] != 0) nv++;
      }
      if ((nv > dfmax) | (tot_iter == max_iter)) {
        for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
        break;
      }

      // Determine eligible set
      if (strcmp(penalty, "lasso")==0) cutoff = 2*lam[l] - lam[l-1];
      if (strcmp(penalty, "MCP")==0) cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lam[l-1]);
      if (strcmp(penalty, "SCAD")==0) cutoff = lam[l] + gamma/(gamma-2)*(lam[l] - lam[l-1]);
      for (int j=0; j<p; j++) if (fabs(z[j]) > (cutoff * alpha * m[j])) e2[j] = 1;
    } else {
      // Determine eligible set
      double lmax = 0;
      for (int j=0; j<p; j++) if (fabs(z[j]) > lmax) lmax = fabs(z[j]);
      if (strcmp(penalty, "lasso")==0) cutoff = 2*lam[l] - lmax;
      if (strcmp(penalty, "MCP")==0) cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lmax);
      if (strcmp(penalty, "SCAD")==0) cutoff = lam[l] + gamma/(gamma-2)*(lam[l] - lmax);
      for (int j=0; j<p; j++) if (fabs(z[j]) > (cutoff * alpha * m[j])) e2[j] = 1;
    }

    while (tot_iter < max_iter) {
      while (tot_iter < max_iter) {
        while (tot_iter < max_iter) {
          // Solve over the active set
          INTEGER(iter)[l]++;
          tot_iter++;
          double maxChange = 0;
          for (int j=0; j<p; j++) {
            if (e1[j]) {
              z[j] = crossprod(X, r, n, j)/n + a[j];

              // Update beta_j
              l1 = lam[l] * m[j] * alpha;
              l2 = lam[l] * m[j] * (1-alpha);
              if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(z[j], l1, l2, gamma, 1);
              if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(z[j], l1, l2, gamma, 1);
              if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(z[j], l1, l2, 1);

              // Update r
              double shift = b[l*p+j] - a[j];
              if (shift !=0) {
                for (int i=0;i<n;i++) r[i] -= shift*X[j*n+i];
                if (fabs(shift) > maxChange) maxChange = fabs(shift);
              }
            }
          }

          // Check for convergence
          for (int j=0; j<p; j++) a[j] = b[l*p+j];
          if (maxChange < eps*sdy) break;
        }

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
    REAL(loss)[l] = g_loss(r, n);
    for (int i=0; i<n; i++) REAL(Eta)[n*l+i] = y[i] - r[i];
  }
  res = cleanup_gaussian(a, r, e1, e2, z, beta, loss, Eta, iter);
  UNPROTECT(4);
  return(res);
}
