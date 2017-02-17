#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
double crossprod(double *X, double *y, int n, int j);
<<<<<<< HEAD
int checkConvergence(double *beta, double *beta_old, double eps, int l, int J);
=======
double wcrossprod(double *X, double *y, double *w, int n, int j);
double S(double z, double l);
>>>>>>> 77281810ac7807a11de80ecd0c0bb5ad70dd09a8
double MCP(double z, double l1, double l2, double gamma, double v);
double SCAD(double z, double l1, double l2, double gamma, double v);
double lasso(double z, double l1, double l2, double v);
double gLoss(double *r, int n);
double sqsum(double *X, int n, int j);

// Memory handling, output formatting (raw)
SEXP cleanupR(double *a, double *r, double *v, double *z, int *e1, int *e2, SEXP beta0, SEXP beta, SEXP loss, SEXP iter) {
  Free(a);
  Free(r);
  Free(v);
  Free(z);
  Free(e1);
  Free(e2);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 4));
  SET_VECTOR_ELT(res, 0, beta0);
  SET_VECTOR_ELT(res, 1, beta);
  SET_VECTOR_ELT(res, 2, loss);
  SET_VECTOR_ELT(res, 3, iter);
  UNPROTECT(5);
  return(res);
}

// Coordinate descent for gaussian models
SEXP cdfit_raw(SEXP X_, SEXP y_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_, SEXP dfmax_, SEXP user_) {

  // Declarations
  int n = length(y_);
  int p = length(X_)/n;
  int L = length(lambda);
  SEXP res, beta0, beta, loss, iter;
  PROTECT(beta0 = allocVector(REALSXP, L));
  double *b0 = REAL(beta0);
  for (int i=0; i<L; i++) b0[i] = 0;
  PROTECT(beta = allocVector(REALSXP, L*p));
  double *b = REAL(beta);
  for (int j=0; j<(L*p); j++) b[j] = 0;
  PROTECT(loss = allocVector(REALSXP, L));
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;
  double *a = Calloc(p, double); // Beta from previous iteration
  for (int j=0; j<p; j++) a[j]=0;
  double a0 = 0; // Beta0 from previous iteration, initially 0 from KKT since y is mean-centered
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
  double *v = Calloc(p, double);
  for (int j=0; j<p; j++) v[j] = sqsum(X, n, j)/n;
<<<<<<< HEAD
  double *z = Calloc(p, double);
  for (int j=0; j<p; j++) z[j] = crossprod(X, r, n, j)/n; // initial a[j] = 0
  int *e1 = Calloc(p, int);
  for (int j=0; j<p; j++) e1[j] = 0;
  int *e2 = Calloc(p, int);
  for (int j=0; j<p; j++) e2[j] = 0;
  double cutoff, l1, l2, mean_resid, shift;
  int converged, lstart;
=======
  int *e = Calloc(p, int);
  for (int j=0; j<p; j++) e[j] = 0;
  double l1, l2, u;
  int lstart;
>>>>>>> 77281810ac7807a11de80ecd0c0bb5ad70dd09a8

  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  double rss = gLoss(r,n);
  if (user) {
    lstart = 0;
  } else {
    REAL(loss)[0] = gLoss(r,n);
    lstart = 1;
  }
  double sdy = sqrt(rss/n);

  // Path
  for (int l=lstart;l<L;l++) {
    R_CheckUserInterrupt();
    if (l != 0) {
      // Assign a0, a
      a0 = b0[l-1];
      for (int j=0;j<p;j++) a[j] = b[(l-1)*p+j];

      // Check dfmax
      int nv = 0;
      for (int j=0; j<p; j++) {
        if (a[j] != 0) nv++;
      }
      if (nv > dfmax) {
      	for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
      	res = cleanupR(a, r, v, z, e1, e2, beta0, beta, loss, iter);
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
<<<<<<< HEAD
      	while (INTEGER(iter)[l] < max_iter) {
      	  // Solve over the active set
      	  INTEGER(iter)[l]++;
      	  // intercept
        	mean_resid = 0.0;
        	for (int i=0; i<n; i++) mean_resid += r[i];
        	mean_resid /= n;
        	b0[l] = mean_resid + a0;
        	for (int i=0; i<n; i++) r[i] -= mean_resid;
      	
      	  for (int j=0; j<p; j++) {
      	    if (e1[j]) {
      	      z[j] = crossprod(X, r, n, j)/n + v[j]*a[j];
      
      	      // Update beta_j
      	      l1 = lam[l] * m[j] * alpha;
      	      l2 = lam[l] * m[j] * (1-alpha);
      	      if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(z[j], l1, l2, gamma, v[j]);
      	      if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(z[j], l1, l2, gamma, v[j]);
      	      if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(z[j], l1, l2, v[j]);
      
      	      // Update r
      	      shift = b[l*p+j] - a[j];
      	      if (shift !=0) for (int i=0;i<n;i++) r[i] -= shift*X[j*n+i];
      	    }
      	  }
=======
	INTEGER(iter)[l]++;
        double maxChange = 0;
	for (int j=0; j<p; j++) {
	  if (e[j]) {
	    u = crossprod(X, r, n, j)/n + v[j]*a[j];
>>>>>>> 77281810ac7807a11de80ecd0c0bb5ad70dd09a8

      	  // Check for convergence
      	  converged = checkConvergence(b, a, eps, l, p);
        	a0 = b0[l];
      	  for (int j=0; j<p; j++) a[j] = b[l*p+j];
      	  if (converged) break;
      	}

<<<<<<< HEAD
      	// Scan for violations in strong set
      	int violations = 0;
      	for (int j=0; j<p; j++) {
      	  if (e1[j]==0 && e2[j]==1) {
      
      	    z[j] = crossprod(X, r, n, j)/n; // a[j] = 0

      	    // Update beta_j
      	    l1 = lam[l] * m[j] * alpha;
      	    l2 = lam[l] * m[j] * (1-alpha);
      	    if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(z[j], l1, l2, gamma, v[j]);
      	    if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(z[j], l1, l2, gamma, v[j]);
      	    if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(z[j], l1, l2, v[j]);
      
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
=======
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
>>>>>>> 77281810ac7807a11de80ecd0c0bb5ad70dd09a8
      }

      // Scan for violations in rest
      int violations = 0;
      for (int j=0; j<p; j++) {
      	if (e2[j]==0) {
      
      	  z[j] = crossprod(X, r, n, j)/n; // a[j] = 0
      
      	  // Update beta_j
      	  l1 = lam[l] * m[j] * alpha;
      	  l2 = lam[l] * m[j] * (1-alpha);
      	  if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(z[j], l1, l2, gamma, v[j]);
      	  if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(z[j], l1, l2, gamma, v[j]);
      	  if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(z[j], l1, l2, v[j]);
      
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
  res = cleanupR(a, r, v, z, e1, e2, beta0, beta, loss, iter);
  return(res);
}
