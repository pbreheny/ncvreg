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

// Memory handling, output formatting (Cox)
SEXP cleanupCox(int *e, double *eta, double *haz, double *rsk, SEXP beta, SEXP Loss, SEXP iter, SEXP residuals, SEXP weights) {
  Free(e);
  Free(eta);
  Free(haz);
  Free(rsk);
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 5));
  SET_VECTOR_ELT(res, 0, beta);
  SET_VECTOR_ELT(res, 1, Loss);
  SET_VECTOR_ELT(res, 2, iter);
  SET_VECTOR_ELT(res, 3, residuals);
  SET_VECTOR_ELT(res, 4, weights);
  UNPROTECT(6);
  return(res);
}

// Coordinate descent for Cox models
SEXP cdfit_cox_dh(SEXP X_, SEXP y_, SEXP d_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_, SEXP dfmax_, SEXP user_, SEXP warn_) {

  // Declarations
  int n = length(y_);
  int p = length(X_)/n;
  int L = length(lambda);
  SEXP res, beta, Loss, iter, residuals, weights;
  PROTECT(beta = allocVector(REALSXP, L*p));
  double *b = REAL(beta);
  for (int j=0; j<(L*p); j++) b[j] = 0;
  PROTECT(residuals = allocVector(REALSXP, n));
  double *r = REAL(residuals);
  PROTECT(weights = allocVector(REALSXP, n));
  double *h = REAL(weights);
  PROTECT(Loss = allocVector(REALSXP, L));
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;
  double *a = Calloc(p, double);    // Beta from previous iteration
  for (int j=0; j<p; j++) a[j] = 0;
  double *X = REAL(X_);
  double *y = REAL(y_);
  double *d = REAL(d_);
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
  double *haz = Calloc(n, double);
  double *rsk = Calloc(n, double);
  int *e = Calloc(p, int);
  for (int j=0; j<p; j++) e[j] = 1;
  double *eta = Calloc(n, double);
  double xwr, xwx, mu, u, v, cutoff, l1, l2, shift, si, exp_eta, s, w;
  int converged, lstart;

  // Initialization
  for (int i=0; i<n; i++) eta[i] = 0;
  // Setup z?  Crossprod at 0?
  // for (int j=0; j<p; j++) z[j] = crossprod(X, s, n, j)/n;

  // If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
  if (user) {
    lstart = 0;
  } else {
    lstart = 1;
    // REAL(Loss)[0] = nullDev; ??
  }

  // Path
  for (int l=lstart; l<L; l++) {
    if (l != 0) {
      // Assign a
      for (int j=0; j<p; j++) a[j] = b[(l-1)*p+j];

      // Check dfmax
      int nv = 0;
      for (int j=0; j<p; j++) {
	if (a[j] != 0) nv++;
      }
      if (nv > dfmax) {
	for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER;
	res = cleanupCox(e, eta, haz, rsk, beta, Loss, iter, residuals, weights);
	return(res);
      }

      // Determine eligible set
      //if (strcmp(penalty, "lasso")==0) cutoff = 2*lam[l] - lam[l-1];
      //if (strcmp(penalty, "MCP")==0) cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lam[l-1]);
      //if (strcmp(penalty, "SCAD")==0) cutoff = lam[l] + gamma/(gamma-2)*(lam[l] - lam[l-1]);
      //for (int j=0; j<p; j++) if (fabs(z[j]) > (cutoff * alpha * m[j])) e2[j] = 1;
    }// else {

      // Determine eligible set
      /* double lmax = 0; */
      /* for (int j=0; j<p; j++) if (fabs(z[j]) > lmax) lmax = fabs(z[j]); */
      /* if (strcmp(penalty, "lasso")==0) cutoff = 2*lam[l] - lmax; */
      /* if (strcmp(penalty, "MCP")==0) cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lmax); */
      /* if (strcmp(penalty, "SCAD")==0) cutoff = lam[l] + gamma/(gamma-2)*(lam[l] - lmax); */
      /* for (int j=0; j<p; j++) if (fabs(z[j]) > (cutoff * alpha * m[j])) e2[j] = 1; */
    //    }

    while (INTEGER(iter)[l] < max_iter) {
      while (INTEGER(iter)[l] < max_iter) {
	while (INTEGER(iter)[l] < max_iter) {
	  INTEGER(iter)[l]++;
	  REAL(Loss)[l] = 0;
	  // Calculate haz, risk
	  for (int i=0; i<n; i++) haz[i] = exp(eta[i]);
	  rsk[n-1] = haz[n-1];
	  for (int i=n-2; i>=0; i--) {
	    rsk[i] = rsk[i+1] + haz[i];
	  }
	  //Rprintf("haz[1]=%f, haz[2]=%f, haz[3]=%f\n", haz[0], haz[1], haz[2]);
	  //Rprintf("rsk[1]=%f, rsk[2]=%f, rsk[3]=%f\n", rsk[0], rsk[1], rsk[2]);
	  // Calculate h, r
	  for (int j=0; j<n; j++) {
	    h[j] = 0;
	    s = d[j];
	    for (int i=0; i<=j; i++) {
	      w = haz[j]/rsk[i];
	      //Rprintf("w: %f, ", w);
	      h[j] += d[i]*w*(1-w);
	      s -= d[i]*w;
	    }
	    //Rprintf("\n");
	    if (h[j]==0) r[j]=0;
	    else r[j] = s/h[j];
	  }
	  //Rprintf("h[1]=%f, h[2]=%f, h[3]=%f\n", h[0], h[1], h[2]);
	  //Rprintf("r[1]=%f, r[2]=%f, r[3]=%f\n", r[0], r[1], r[2]);
	  /*   mu = exp(eta[i]); */
	  /*   w[i] = mu; */
	  /*   s[i] = y[i] - mu; */
	  /*   r[i] = s[i]/w[i]; */
	  /*   if (y[i]!=0) REAL(Dev)[l] += y[i]*log(y[i]/mu); */
	  /* } */
	  /* if (REAL(Dev)[l]/nullDev < .01) { */
	  /*   if (warn) warning("Model saturated; exiting..."); */
	  /*   for (int ll=l; ll<L; ll++) INTEGER(iter)[ll] = NA_INTEGER; */
	  /*   res = cleanupP(s, w, a, r, e1, e2, z, eta, beta0, beta, Dev, iter); */
	  /*   return(res); */
	  /* } */

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
	      //Rprintf("u=%f, v=%f, b=%f\n", u, v, b[l*p+j]);

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
	  for (int j=0; j<p; j++) a[j] = b[l*p+j];
	  if (converged) break;
	}

	// Scan for violations in strong set
	int violations = 0;
	/* for (int j=0; j<p; j++) { */
	/*   if (e1[j]==0 & e2[j]==1) { */
	/*     z[j] = crossprod(X, s, n, j)/n; */
	/*     l1 = lam[l] * m[j] * alpha; */
	/*     if (fabs(z[j]) > l1) { */
	/*       e1[j] = e2[j] = 1; */
	/*       violations++; */
	/*     } */
	/*   } */
	/* } */
	if (violations==0) break;
      }

      // Scan for violations in rest
      int violations = 0;
      /* for (int j=0; j<p; j++) { */
      /* 	if (e2[j]==0) { */
      /* 	  z[j] = crossprod(X, s, n, j)/n; */
      /* 	  l1 = lam[l] * m[j] * alpha; */
      /* 	  if (fabs(z[j]) > l1) { */
      /* 	    e1[j] = e2[j] = 1; */
      /* 	    violations++; */
      /* 	  } */
      /* 	} */
      /* } */
      if (violations==0) break;
    }
  }
  res = cleanupCox(e, eta, haz, rsk, beta, Loss, iter, residuals, weights);
  return(res);
}
