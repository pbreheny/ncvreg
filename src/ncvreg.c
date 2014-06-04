#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>
SEXP cdfit_gaussian(SEXP X_, SEXP y_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_, SEXP dfmax_, SEXP user_);
SEXP cdfit_binomial(SEXP X_, SEXP y_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_, SEXP dfmax_, SEXP user_, SEXP warn_);
SEXP cdfit_poisson(SEXP X_, SEXP y_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_, SEXP dfmax_, SEXP user_, SEXP warn_);
SEXP cdfit_cox_dh(SEXP X_, SEXP y_, SEXP d_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_, SEXP dfmax_, SEXP user_, SEXP warn_);
SEXP cdfit_raw(SEXP X_, SEXP y_, SEXP penalty_, SEXP lambda, SEXP eps_, SEXP max_iter_, SEXP gamma_, SEXP multiplier, SEXP alpha_, SEXP dfmax_);
SEXP standardize(SEXP X_);
SEXP maxprod(SEXP X_, SEXP y_, SEXP v_, SEXP m_);


// Cross product of y with jth column of X
double crossprod(double *X, double *y, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += X[nn+i]*y[i];
  return(val);
}

// Weighted cross product of y with jth column of x
double wcrossprod(double *X, double *y, double *w, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += X[nn+i]*y[i]*w[i];
  return(val);
}

// Weighted sum of squares of jth column of X
double wsqsum(double *X, double *w, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += w[i] * pow(X[nn+i], 2);
  return(val);
}

// Sum of squares of jth column of X
double sqsum(double *X, int n, int j) {
  int nn = n*j;
  double val=0;
  for (int i=0;i<n;i++) val += pow(X[nn+i], 2);
  return(val);
}

double sum(double *x, int n) {
  double val=0;
  for (int i=0;i<n;i++) val += x[i];
  return(val);
}

int checkConvergence(double *beta, double *beta_old, double eps, int l, int J) {
  int converged = 1;
  for (int j=0; j<J; j++) {
    if (fabs((beta[l*J+j]-beta_old[j])/beta_old[j]) > eps) {
      converged = 0;
      break;
    }
  }
  return(converged);
}

double MCP(double z, double l1, double l2, double gamma, double v) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-l1)/(v*(1+l2-1/gamma)));
  else return(z/(v*(1+l2)));
}

double SCAD(double z, double l1, double l2, double gamma, double v) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else if (fabs(z) <= (l1*(1+l2)+l1)) return(s*(fabs(z)-l1)/(v*(1+l2)));
  else if (fabs(z) <= gamma*l1*(1+l2)) return(s*(fabs(z)-gamma*l1/(gamma-1))/(v*(1-1/(gamma-1)+l2)));
  else return(z/(v*(1+l2)));
}

double lasso(double z, double l1, double l2, double v) {
  double s=0;
  if (z > 0) s = 1;
  else if (z < 0) s = -1;
  if (fabs(z) <= l1) return(0);
  else return(s*(fabs(z)-l1)/(v*(1+l2)));
}

static R_CallMethodDef callMethods[] = {
  {"cdfit_gaussian", (DL_FUNC) &cdfit_gaussian, 11},
  {"cdfit_raw", (DL_FUNC) &cdfit_raw, 10},
  {"cdfit_binomial", (DL_FUNC) &cdfit_binomial, 12},
  {"cdfit_poisson", (DL_FUNC) &cdfit_poisson, 12},
  {"cdfit_cox_dh", (DL_FUNC) &cdfit_cox_dh, 13},
  {"standardize", (DL_FUNC) &standardize, 1},
  {"maxprod", (DL_FUNC) &maxprod, 4},
  {NULL, NULL, 0}
};

void R_init_ncvreg(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
