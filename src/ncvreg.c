#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>

// Cross product of y with jth column of X
static double crossprod(double *X, double *y, int n, int j) {
  int nn = n*j;
  double val;
  for (int i=0;i<n;i++) val += X[nn+i]*y[i];
  return(val);
}

// Weighted cross product of y with jth column of x
static double wcrossprod(double *X, double *y, double *w, int n, int j) {
  int nn = n*j;
  double val;
  for (int i=0;i<n;i++) val += X[nn+i]*y[i]*w[i];
  return(val);
}

// Weighted sum of squares of jth column of X
static double wsqsum(double *X, double *w, int n, int j) {
  int nn = n*j;
  double val;
  for (int i=0;i<n;i++) val += w[i] * pow(X[nn+i], 2);
  return(val);
}

static double sum(double *x, int n) {
  double val;
  for (int i=0;i<n;i++) val += x[i];
  return(val);
}

static int checkConvergence(double *beta, double *beta_old, double eps, int l, int J, int *active, int nactive) {
  int converged = 1;
  for (int jj=0;jj<nactive;jj++) {
    int j = active[jj];
    if (fabs((beta[l*J+j]-beta_old[j])/beta_old[j]) > eps) {
      converged = 0;
      break;
    }
  }
  return(converged);
}

// Gaussian loss
static double gLoss(double *r, int n)
{
  double l = 0;
  for (int i=0;i<n;i++) l = l + pow(r[i],2);
  return(l);
}

// Soft thresholding
static double S(double z, double l)
{
  if (z > l) return(z-l);
  if (z < -l) return(z+l);
  return(0);
}

static double MCP(double z, double l1, double l2, double gamma, double v)
{
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-l1)/(v*(1+l2-1/gamma)));
  else return(z/(v*(1+l2)));
}

static double SCAD(double z, double l1, double l2, double gamma, double v)
{
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= (l1*(1+l2)+l1)) return(s*(fabs(z)-l1)/(v*(1+l2)));
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-gamma*l1/(gamma-1))/(v*(1-1/(gamma-1)+l2)));
  else return(z/(v*(1+l2)));
}

static double lasso(double z, double l1, double l2, double v)
{
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else return(s*(fabs(z)-l1)/(v*(1+l2)));
}

// Coordinate descent for gaussian models
SEXP cdfit_gaussian(SEXP X_, SEXP y_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_, SEXP dfmax_, SEXP user_) {

  // Declarations
  int n = length(y_);
  int p = length(X_)/n;
  int L = length(lambda);
  SEXP res, beta, loss, iter;
  PROTECT(res = allocVector(VECSXP, 3));
  PROTECT(beta = allocVector(REALSXP, L*p));
  double *b = REAL(beta);
  for (int i=0; i<(L*p); i++) b[i] = 0;
  PROTECT(loss = allocVector(REALSXP, L));
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;
  double *a = Calloc(p, double); // Beta from previous iteration
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
  int converged, lstart;
  double *r = Calloc(n, double);
  for (int i=0;i<n;i++) r[i] = y[i];
  int *active = Calloc(p, int);
  int nactive = 0;

  if (user) {
    lstart = 0;
  } else {
    REAL(loss)[0] = gLoss(r,n);
    lstart = 1;
  }

  // Path
  for (int l=lstart;l<L;l++) {
    if (l != 0) for (int j=0;j<p;j++) a[j] = b[(l-1)*p+j];
    while (INTEGER(iter)[l] < max_iter) {

      // Check dfmax
      if (nactive > dfmax) {
	for (int ll=l;ll<L;ll++) {
	  for (int j=0;j<p;j++) b[ll*p+j] = R_NaReal;
	}
	Free(a);
	Free(r);
	Free(active);
	SET_VECTOR_ELT(res, 0, beta);
	SET_VECTOR_ELT(res, 1, loss);
	SET_VECTOR_ELT(res, 2, iter);
	UNPROTECT(4);
	return(res);
      }

      // Active set
      if (nactive > 0) {
	while (INTEGER(iter)[l] < max_iter) {
	  INTEGER(iter)[l]++;

	  // Covariates
	  for (int jj=0;jj<nactive;jj++) {
	    int j = active[jj];

	    // Calculate z
	    double z = crossprod(X, r, n, j)/n + a[j];

	    // Update beta_j
	    double l1 = lam[l] * m[j] * alpha;
	    double l2 = lam[l] * m[j] * (1-alpha);
	    if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(z, l1, l2, gamma, 1);
	    if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(z, l1, l2, gamma, 1);
	    if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(z, l1, l2, 1);

	    // If dropout, adjust active
	    if (b[l*p+j]==0) {
	      //Rprintf("DROPOUT\n");
	      for (int jjj=jj;jjj<nactive;jjj++) active[jjj] = active[jjj+1];
	      jj--;
	      nactive--;
	    }

	    // Update r
	    double shift = b[l*p+j] - a[j];
            for (int i=0;i<n;i++) r[i] -= shift*X[j*n+i];
	  }

	  // Check for convergence
	  converged = checkConvergence(b, a, eps, l, p, active, nactive);
	  for (int jj=0;jj<nactive;jj++) {
	    int j = active[jj];
	    a[j] = b[l*p+j];
	  }
	  if (converged) break;
	}
      }

      // Scan for new variables
      int newvar = 0;
      for (int j=0; j<p; j++) {
	if (b[l*p+j] == 0) {
	  double z = crossprod(X, r, n, j)/n + a[j];

	  // Update beta_j
	  double l1 = lam[l] * m[j] * alpha;
	  double l2 = lam[l] * m[j] * (1-alpha);
	  if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(z, l1, l2, gamma, 1);
	  if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(z, l1, l2, gamma, 1);
	  if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(z, l1, l2, 1);

	  // If something enters the active set, update active set, residuals
	  if (b[l*p+j] !=0) {
	    active[nactive] = j;
	    nactive++;
	    double shift = b[l*p+j] - a[j];
	    for (int i=0;i<n;i++) r[i] -= shift*X[j*n+i];
	    a[j] = b[l*p+j];
	    newvar=1;
	  }
	}
      }
      //Rprintf("nactive: %d\n", nactive);
      //Rprintf("active[0]: %f\n", active[0]);

      if (newvar==0) {
	REAL(loss)[l] = gLoss(r,n);
	break;
      }
      //for (int j=0;j<p;j++) a[j] = b[l*p+j];
    }
  }
  Free(a);
  Free(r);
  Free(active);
  SET_VECTOR_ELT(res, 0, beta);
  SET_VECTOR_ELT(res, 1, loss);
  SET_VECTOR_ELT(res, 2, iter);
  UNPROTECT(4);
  return(res);
}

// Coordinate descent for binomial models
SEXP cdfit_binomial(SEXP X_, SEXP y_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_, SEXP dfmax_, SEXP user_, SEXP warn_) {

  // Declarations
  int n = length(y_);
  int p = length(X_)/n;
  int L = length(lambda);
  SEXP res, beta0, beta, Dev, iter;
  PROTECT(res = allocVector(VECSXP, 4));
  PROTECT(beta0 = allocVector(REALSXP, L));
  double *b0 = REAL(beta0);
  for (int i=0; i<L; i++) b0[i] = 0;
  PROTECT(beta = allocVector(REALSXP, L*p));
  double *b = REAL(beta);
  for (int i=0; i<(L*p); i++) b[i] = 0;
  PROTECT(Dev = allocVector(REALSXP, L));
  PROTECT(iter = allocVector(INTSXP, L));
  for (int i=0; i<L; i++) INTEGER(iter)[i] = 0;
  double *a = Calloc(p, double); // Beta from previous iteration
  double a0;                    // Beta0 from previous iteration
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
  int converged, lstart;
  double *r = Calloc(n, double);
  double *w = Calloc(n, double);
  int *active = Calloc(p, int);
  int nactive = 0;

  // Initialization
  double ybar=0;
  for (int i=0;i<n;i++) ybar = ybar + y[i];
  ybar = ybar/n;
  if (user) a0 = log(ybar/(1-ybar));
  else b0[0] = log(ybar/(1-ybar));
  double nullDev = 0;
  for (int i=0;i<n;i++) nullDev = nullDev - y[i]*log(ybar)-(1-y[i])*log(1-ybar);

  if (user) lstart = 0;
  else {
    lstart = 1;
    REAL(Dev)[0] = nullDev;
  }

  // Path
  double xwr,xwx,eta,pi,z,v;
  for (int l=lstart;l<L;l++) {
    if (l != 0) {
      a0 = b0[l-1];
      for (int j=0;j<p;j++) a[j] = b[(l-1)*p+j];
    }

    while (INTEGER(iter)[l] < max_iter) {
      //Rprintf("nactive: %d\n", nactive);

      // Check dfmax
      if (nactive > dfmax) {
	for (int ll=l;ll<L;ll++) {
	  for (int j=0;j<p;j++) b[ll*p+j] = R_NaReal;
	}
	Free(a);
	Free(w);
	Free(r);
	Free(active);
	SET_VECTOR_ELT(res, 0, beta0);
	SET_VECTOR_ELT(res, 1, beta);
	SET_VECTOR_ELT(res, 2, Dev);
	SET_VECTOR_ELT(res, 3, iter);
	UNPROTECT(5);
	return(res);
      }

      // Active set
      while (INTEGER(iter)[l] < max_iter) {
	INTEGER(iter)[l]++;

	// Approximate L
	REAL(Dev)[l] = 0;
	for (int i=0;i<n;i++) {
	  eta = a0;
	  for (int j=0;j<p;j++) eta = eta + X[j*n+i]*a[j];
	  pi = exp(eta)/(1+exp(eta));
	  if (pi > .9999) {
	    pi = 1;
	    w[i] = .0001;
	  }
	  else if (pi < .0001) {
	    pi = 0;
	    w[i] = .0001;
	  }
	  else w[i] = pi*(1-pi);
	  r[i] = (y[i] - pi)/w[i];
	  if (y[i]==1) REAL(Dev)[l] = REAL(Dev)[l] - log(pi);
	  if (y[i]==0) REAL(Dev)[l] = REAL(Dev)[l] - log(1-pi);
	}
	if (REAL(Dev)[l]/nullDev < .01) {
	  if (warn) warning("Model saturated; exiting...");
	  for (int ll=l;ll<L;ll++) {
	    b0[ll] = R_NaReal;
	    for (int j=0;j<p;j++) b[ll*p+j] = R_NaReal;
	  }
	  Free(a);
	  Free(w);
	  Free(r);
	  Free(active);
	  SET_VECTOR_ELT(res, 0, beta0);
	  SET_VECTOR_ELT(res, 1, beta);
	  SET_VECTOR_ELT(res, 2, Dev);
	  SET_VECTOR_ELT(res, 3, iter);
	  UNPROTECT(5);
	  return(res);
	}

	// Intercept
	xwr = crossprod(w, r, n, 0);
	xwx = sum(w, n);
	b0[l] = xwr/xwx + a0;
	for (int i=0;i<n;i++) r[i] = r[i] - (b0[l] - a0);

	// Covariates
	if (nactive > 0) {
	  for (int jj=0;jj<nactive;jj++) {
	    int j = active[jj];

	    // Calculate z
	    xwr = wcrossprod(X, r, w, n, j);
	    xwx = wsqsum(X, w, n, j);
	    z = xwr/n + (xwx/n)*a[j];
	    v = xwx/n;

	    // Update b_j
	    double l1 = lam[l] * m[j] * alpha;
	    double l2 = lam[l] * m[j] * (1-alpha);
	    if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(z, l1, l2, gamma, v);
	    if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(z, l1, l2, gamma, v);
	    if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(z, l1, l2, v);

	    // If dropout, adjust active
	    if (b[l*p+j]==0) {
	      // Rprintf("DROPOUT!!!!\n");
	      for (int jjj=jj;jjj<nactive;jjj++) active[jjj] = active[jjj+1];
	      jj--;
	      nactive--;
	    }

	    // Update r
	    double shift = b[l*p+j] - a[j];
	    for (int i=0;i<n;i++) r[i] -= shift*X[j*n+i];
	  }
	}

	// Check for convergence
	converged = checkConvergence(b, a, eps, l, p, active, nactive);
	a0 = b0[l];
	for (int jj=0;jj<nactive;jj++) {
	  int j = active[jj];
	  a[j] = b[l*p+j];
	}
	if (converged) break;
      }

      // Scan for new variables
      int newvar = 0;
      for (int j=0; j<p; j++) {
	if (b[l*p+j] == 0) {
	  // Calculate z
	  xwr = wcrossprod(X, r, w, n, j);
	  xwx = wsqsum(X, w, n, j);
	  z = xwr/n + (xwx/n)*a[j];
	  v = xwx/n;

	  // Update b_j
	  double l1 = lam[l] * m[j] * alpha;
	  double l2 = lam[l] * m[j] * (1-alpha);
	  if (strcmp(penalty,"MCP")==0) b[l*p+j] = MCP(z, l1, l2, gamma, v);
	  if (strcmp(penalty,"SCAD")==0) b[l*p+j] = SCAD(z, l1, l2, gamma, v);
	  if (strcmp(penalty,"lasso")==0) b[l*p+j] = lasso(z, l1, l2, v);

	  // If something enters the active set, update active set, residuals
	  if (b[l*p+j] !=0) {
	    active[nactive] = j;
	    nactive++;
	    double shift = b[l*p+j] - a[j];
	    for (int i=0;i<n;i++) r[i] -= shift*X[j*n+i];
	    a[j] = b[l*p+j];
	    newvar=1;
	  }
	}
      }

      if (newvar==0) break;
    }
  }
  Free(a);
  Free(w);
  Free(r);
  Free(active);
  SET_VECTOR_ELT(res, 0, beta0);
  SET_VECTOR_ELT(res, 1, beta);
  SET_VECTOR_ELT(res, 2, Dev);
  SET_VECTOR_ELT(res, 3, iter);
  UNPROTECT(5);
  return(res);
}

/* static const R_CMethodDef cMethods[] = { */
/*   {"cdfit_gaussian", (DL_FUNC) &cdfit_gaussian, 17}, */
/*   {"cdfit_binomial", (DL_FUNC) &cdfit_binomial, 19}, */
/*   NULL */
/* }; */

static R_CallMethodDef callMethods[] = {
  {"cdfit_gaussian", (DL_FUNC) &cdfit_gaussian, 11},
  {"cdfit_binomial", (DL_FUNC) &cdfit_binomial, 12},
  NULL
};

void R_init_ncvreg(DllInfo *info)
{
  R_registerRoutines(info,NULL,callMethods,NULL,NULL);
}
