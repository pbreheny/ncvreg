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
int checkConvergence(double *beta, double *beta_old, double eps, int l, int J);
double S(double z, double l);
double MCP(double z, double l1, double l2, double gamma, double v);
double SCAD(double z, double l1, double l2, double gamma, double v);
double lasso(double z, double l1, double l2, double v);

// Memory handling, output formatting (Gaussian)
SEXP cleanupG(double *a, double *r, int *e, double *z, SEXP beta, SEXP loss, SEXP iter) {
  Free(a);
  Free(r);
  Free(e);
  Free(z);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(res, 0, beta);
  SET_VECTOR_ELT(res, 1, loss);
  SET_VECTOR_ELT(res, 2, iter);
  UNPROTECT(4);
  return(res);
}

// Gaussian loss
double gLoss(double *r, int n)
{
  double l = 0;
  for (int i=0;i<n;i++) l = l + pow(r[i],2);
  return(l);
}

// Coordinate descent for gaussian models
SEXP cdfit_gaussian(SEXP X_, SEXP y_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_, SEXP dfmax_, SEXP user_) {

  // Declarations
  int n = length(y_);
  int p = length(X_)/n;
  int L = length(lambda);
  SEXP res, beta, loss, iter;
  PROTECT(beta = allocVector(REALSXP, L*p));
  double *b = REAL(beta);
  for (int j=0; j<(L*p); j++) b[j] = 0;
  PROTECT(loss = allocVector(REALSXP, L));
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;
  double *a = Calloc(p, double); // Beta from previous iteration
  for (int j=0; j<p; j++) a[j]=0;
  double *X = REAL(X_);
  double *y = REAL(y_);
  const char *penalty = CHAR(STRING_ELT(penalty_, 0));
  double *lam = REAL(lambda);
  double eps = REAL(eps_)[0];
  int max_iter = INTEGER(max_iter_)[0];
  double gamma = REAL(gamma_)[0];
  double *m = REAL(multiplier);
  double alpha = REAL(alpha_)[0];
  int dfmax = INTEGER(dfmax_)[0];
  int user = INTEGER(user_)[0];
  double *r = Calloc(n, double);
  for (int i=0; i<n; i++) r[i] = y[i];
  double *z = Calloc(p, double);
  for (int j=0; j<p; j++) z[j] = crossprod(X, r, n, j)/n;
  int *e = Calloc(p, int);
  for (int j=0; j<p; j++) e[j] = 0;
  double cutoff, l1, l2;
  int converged, lstart;

  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  if (user) {
    lstart = 0;
  } else {
    REAL(loss)[0] = gLoss(r,n);
    lstart = 1;
  }

  // Path
  for (int l=lstart;l<L;l++) {
    if (l != 0) {
      // Assign a
      for (int j=0;j<p;j++) a[j] = b[(l-1)*p+j];

      // Check dfmax
      int nv = 0;
      for (int j=0; j<p; j++) {
	if (a[j] != 0) nv++;
      }
      if (nv > dfmax) {
	for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
	res = cleanupG(a, r, e, z, beta, loss, iter);
	return(res);
      }

      // Determine eligible set
      if (strcmp(penalty, "lasso")==0) cutoff = 2*lam[l] - lam[l-1];
      if (strcmp(penalty, "MCP")==0) cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lam[l-1]);
      if (strcmp(penalty, "SCAD")==0) cutoff = lam[l] + gamma/(gamma-2)*(lam[l] - lam[l-1]);
      for (int j=0; j<p; j++) if (fabs(z[j]) > (cutoff * alpha * m[j])) e[j] = 1;
    } else {
      // Determine eligible set
      double lmax = 0;
      for (int j=0; j<p; j++) if (fabs(z[j]) > lmax) lmax = fabs(z[j]);
      if (strcmp(penalty, "lasso")==0) cutoff = 2*lam[l] - lmax;
      if (strcmp(penalty, "MCP")==0) cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lmax);
      if (strcmp(penalty, "SCAD")==0) cutoff = lam[l] + gamma/(gamma-2)*(lam[l] - lmax);
      for (int j=0; j<p; j++) if (fabs(z[j]) > (cutoff * alpha * m[j])) e[j] = 1;
    }

    while (INTEGER(iter)[l] < max_iter) {

      // Solve over eligible set
      while (INTEGER(iter)[l] < max_iter) {
	INTEGER(iter)[l]++;
	for (int j=0; j<p; j++) {
	  if (e[j]) {
	    // Rprintf("a: %f\n", a[j]);
	    z[j] = crossprod(X, r, n, j)/n + a[j];

	    // Update beta_j
	    l1 = lam[l] * m[j] * alpha;
	    l2 = lam[l] * m[j] * (1-alpha);
	    if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(z[j], l1, l2, gamma, 1);
	    if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(z[j], l1, l2, gamma, 1);
	    if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(z[j], l1, l2, 1);
	    // if (j==0) Rprintf("b[0] = %f\n", b[l*p+j]);

	    // Update r
	    double shift = b[l*p+j] - a[j];
	    if (shift !=0) for (int i=0;i<n;i++) r[i] -= shift*X[j*n+i];
	  }
	}

	// Check for convergence
	converged = checkConvergence(b, a, eps, l, p);
	for (int j=0; j<p; j++) a[j] = b[l*p+j];
	if (converged) break;
      }

      // Scan for violations
      int violations = 0;
      for (int j=0; j<p; j++) {
	if (e[j]==0) {

	  z[j] = crossprod(X, r, n, j)/n;

	  // Update beta_j
	  l1 = lam[l] * m[j] * alpha;
	  l2 = lam[l] * m[j] * (1-alpha);
	  if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(z[j], l1, l2, gamma, 1);
	  if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(z[j], l1, l2, gamma, 1);
	  if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(z[j], l1, l2, 1);

	  // If something enters the eligible set, update eligible set & residuals
	  if (b[l*p+j] !=0) {
	    e[j] = 1;
	    for (int i=0; i<n; i++) r[i] -= b[l*p+j]*X[j*n+i];
	    a[j] = b[l*p+j];
	    violations++;
	  }
	}
      }

      if (violations==0) {
	REAL(loss)[l] = gLoss(r, n);
	break;
      }
    }
  }
  res = cleanupG(a, r, e, z, beta, loss, iter);
  return(res);
}
