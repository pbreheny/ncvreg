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
double MCP(double z, double l1, double l2, double gamma, double v);
double SCAD(double z, double l1, double l2, double gamma, double v);
double lasso(double z, double l1, double l2, double v);

// Memory handling, output formatting (Poisson)
SEXP cleanupP(double *s, double *w, double *a, double *r, int *e1, int *e2, double *z, double *eta, SEXP beta0, SEXP beta, SEXP Dev, SEXP iter) {
  Free(s);
  Free(w);
  Free(a);
  Free(r);
  Free(e1);
  Free(e2);
  Free(z);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(res, 0, beta0);
  SET_VECTOR_ELT(res, 1, beta);
  SET_VECTOR_ELT(res, 2, Dev);
  SET_VECTOR_ELT(res, 3, iter);
  UNPROTECT(5);
  return(res);
}

// Coordinate descent for poisson models
SEXP cdfit_poisson(SEXP X_, SEXP y_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_, SEXP dfmax_, SEXP user_, SEXP warn_) {

  // Declarations
  int n = length(y_);
  int p = length(X_)/n;
  int L = length(lambda);
  SEXP res, beta0, beta, Dev, iter;
  PROTECT(beta0 = allocVector(REALSXP, L));
  double *b0 = REAL(beta0);
  for (int i=0; i<L; i++) b0[i] = 0;
  PROTECT(beta = allocVector(REALSXP, L*p));
  double *b = REAL(beta);
  for (int j=0; j<(L*p); j++) b[j] = 0;
  PROTECT(Dev = allocVector(REALSXP, L));
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;
  double *a = Calloc(p, double);    // Beta from previous iteration
  for (int j=0; j<p; j++) a[j] = 0;
  double a0 = 0;                    // Beta0 from previous iteration
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
  int warn = INTEGER(warn_)[0];
  double *r = Calloc(n, double);
  double *w = Calloc(n, double);
  double *s = Calloc(n, double);
  double *z = Calloc(p, double);
  double *eta = Calloc(n, double);
  int *e1 = Calloc(p, int);
  for (int j=0; j<p; j++) e1[j] = 0;
  int *e2 = Calloc(p, int);
  for (int j=0; j<p; j++) e2[j] = 0;
  double xwr, xwx, mu, u, v, cutoff, l1, l2, shift, si;
  int converged, lstart;

  // Initialization
  double ybar = sum(y, n)/n;
  a0 = b0[0] = log(ybar);
  double nullDev = 0;
  for (int i=0;i<n;i++) if (y[i]!=0) nullDev += y[i]*log(y[i]/ybar);
  for (int i=0; i<n; i++) s[i] = y[i] - ybar;
  for (int i=0; i<n; i++) eta[i] = a0;
  for (int j=0; j<p; j++) z[j] = crossprod(X, s, n, j)/n;

  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  if (user) {
    lstart = 0;
  } else {
    lstart = 1;
    REAL(Dev)[0] = nullDev;
  }

  // Path
  for (int l=lstart; l<L; l++) {
    if (l != 0) {
      // Assign a, a0
      a0 = b0[l-1];
      for (int j=0; j<p; j++) a[j] = b[(l-1)*p+j];

      // Check dfmax
      int nv = 0;
      for (int j=0; j<p; j++) {
	if (a[j] != 0) nv++;
      }
      if (nv > dfmax) {
	for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
	res = cleanupP(s, w, a, r, e1, e2, z, eta, beta0, beta, Dev, iter);
	return(res);
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

    while (INTEGER(iter)[l] < max_iter) {
      while (INTEGER(iter)[l] < max_iter) {
	while (INTEGER(iter)[l] < max_iter) {
	  INTEGER(iter)[l]++;
	  REAL(Dev)[l] = 0;
	  for (int i=0;i<n;i++) {
	    mu = exp(eta[i]);
	    w[i] = mu;
	    s[i] = y[i] - mu;
	    r[i] = s[i]/w[i];
	    if (y[i]!=0) REAL(Dev)[l] += y[i]*log(y[i]/mu);
	  }
	  if (REAL(Dev)[l]/nullDev < .01) {
	    if (warn) warning("Model saturated; exiting...");
	    for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
	    res = cleanupP(s, w, a, r, e1, e2, z, eta, beta0, beta, Dev, iter);
	    return(res);
	  }

	  // Intercept
	  xwr = crossprod(w, r, n, 0);
	  xwx = sum(w, n);
	  b0[l] = xwr/xwx + a0;
	  for (int i=0; i<n; i++) {
	    si = b0[l] - a0;
	    r[i] -= si;
	    eta[i] += si;
	  }

	  // Covariates
	  for (int j=0; j<p; j++) {
	    if (e1[j]) {

	      // Calculate u, v
	      xwr = wcrossprod(X, r, w, n, j);
	      xwx = wsqsum(X, w, n, j);
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
		/* for (int i=0;i<n;i++) r[i] -= shift*X[j*n+i]; */
		/* for (int i=0;i<n;i++) eta[i] += shift*X[j*n+i]; */
		for (int i=0;i<n;i++) {
		  si = shift*X[j*n+i];
		  r[i] -= si;
		  eta[i] += si;
		}
	      }
	    }
	  }

	  // Check for convergence
	  converged = checkConvergence(b, a, eps, l, p);
	  a0 = b0[l];
	  for (int j=0; j<p; j++) a[j] = b[l*p+j];
	  if (converged) break;
	}

	// Scan for violations in strong set
	int violations = 0;
	for (int j=0; j<p; j++) {
	  if (e1[j]==0 & e2[j]==1) {
	    z[j] = crossprod(X, s, n, j)/n;
	    l1 = lam[l] * m[j] * alpha;
	    if (fabs(z[j]) > l1) {
	      e1[j] = e2[j] = 1;
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
	  z[j] = crossprod(X, s, n, j)/n;
	  l1 = lam[l] * m[j] * alpha;
	  if (fabs(z[j]) > l1) {
	    e1[j] = e2[j] = 1;
	    violations++;
	  }
	}
      }
      if (violations==0) break;
    }
  }
  res = cleanupP(s, w, a, r, e1, e2, z, eta, beta0, beta, Dev, iter);
  return(res);
}
