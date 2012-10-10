#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>

static double *vector(int n)
{
  double *v;
  v = Calloc(n, double);
  return v;
}

static void free_vector(double *v)
{
  Free(v);
}

static int checkConvergence(double *beta, double *beta_old, double eps, int l, int J)
{
  int j;
  int converged = 1;
  for (j=0; j < J; j++) {
    if (beta[l*J+j]!=0 & beta_old[j]!=0) {
      if (fabs((beta[l*J+j]-beta_old[j])/beta_old[j]) > eps) {
	converged = 0;
	break;
      }
    } else if (beta[l*J+j]==0 & beta_old[j]!=0) {
      converged = 0;
      break;
    } else if (beta[l*J+j]!=0 & beta_old[j]==0) {
      converged = 0;
      break;
    }
  }
  return(converged);
}

static double gLoss(double *r, int n)
{
  double l = 0;
  for (int i=0;i<n;i++) l = l + pow(r[i],2);
  return(l);
}

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

static void cdfit_gaussian(double *beta, double *loss, int *iter, double *x, double *y, int *n_, int *p_, char **penalty_, double *lambda, int *L_, double *eps_, int *max_iter_, double *gamma_, double *multiplier, double *alpha_, int *dfmax_, int *user_)
{
  /* Declarations */
  int L=L_[0]; int p=p_[0]; int n=n_[0]; int max_iter=max_iter_[0]; double eps=eps_[0]; double gamma=gamma_[0]; double alpha=alpha_[0]; char *penalty=penalty_[0]; int dfmax=dfmax_[0]; int user=user_[0];
  int converged, active, lstart;
  double *r, *beta_old;
  r = vector(n); for (int i=0;i<n;i++) r[i] = y[i];
  beta_old = vector(p);

  if (user) lstart = 0;
  else
    {
      loss[0] = gLoss(r,n);
      lstart = 1;
    }

  /* Path */
  for (int l=lstart;l<L;l++) {
    if (l != 0) for (int j=0;j<p;j++) beta_old[j] = beta[(l-1)*p+j];
    while (iter[l] < max_iter) {
      converged = 0;
      iter[l] = iter[l] + 1;

      /* Check dfmax */
      active = 0;
      for (int j=0;j<p;j++) if (beta[l*p+j]!=0) active++;
      if (active > dfmax) {
	for (int ll=l;ll<L;ll++) {
	  for (int j=0;j<p;j++) beta[ll*p+j] = R_NaReal;
	}
	free_vector(beta_old);
	free_vector(r);
	return;
      }

      /* Covariates */
      for (int j=0;j<p;j++)
	{
	  /* Calculate z */
	  double z = 0;
	  for (int i=0;i<n;i++) z = z + x[j*n+i]*r[i];
	  z = z/n + beta_old[j];

	  /* Update beta_j */
	  double l1 = lambda[l] * multiplier[j] * alpha;
	  double l2 = lambda[l] * multiplier[j] * (1-alpha);
	  if (strcmp(penalty,"MCP")==0) beta[l*p+j] = MCP(z, l1, l2, gamma, 1);
	  if (strcmp(penalty,"SCAD")==0) beta[l*p+j] = SCAD(z, l1, l2, gamma, 1);
	  if (strcmp(penalty,"lasso")==0) beta[l*p+j] = SCAD(z, l1, l2, gamma, 1);

	  /* Update r */
	  if (beta[l*p+j] != beta_old[j]) for (int i=0;i<n;i++) r[i] = r[i] - (beta[l*p+j] - beta_old[j])*x[j*n+i];
	}

      /* Check for convergence */
      if (checkConvergence(beta,beta_old,eps,l,p))
	{
	  converged  = 1;
	  loss[l] = gLoss(r,n);
	  break;
	}
      for (int j=0;j<p;j++) beta_old[j] = beta[l*p+j];
    }
    /*if (converged==0) warning("Failed to converge");*/
  }

  free_vector(beta_old);
  free_vector(r);
}

static void cdfit_binomial(double *beta0, double *beta, double *Dev, int *iter, double *x, double *y, int *n_, int *p_, char **penalty_, double *lambda, int *L_, double *eps_, int *max_iter_, double *gamma_, double *multiplier, double *alpha_, int *dfmax_, int *user_, int *warn_)
{
  /* Declarations */
  int L=L_[0];int p=p_[0];int n=n_[0];int max_iter=max_iter_[0];double eps=eps_[0];double gamma=gamma_[0]; double alpha=alpha_[0];char *penalty=penalty_[0];int dfmax=dfmax_[0]; int user=user_[0]; int warn = warn_[0];
  int converged, active, lstart;
  double beta0_old;
  double *r, *w, *beta_old;
  r = vector(n);
  w = vector(n);
  beta_old = vector(p);

  /* Initialization */
  double ybar=0;
  for (int i=0;i<n;i++) ybar = ybar + y[i];
  ybar = ybar/n;
  if (user) beta0_old = log(ybar/(1-ybar));
  else beta0[0] = log(ybar/(1-ybar));
  double nullDev = 0;
  for (int i=0;i<n;i++) nullDev = nullDev - y[i]*log(ybar)-(1-y[i])*log(1-ybar);

  if (user) lstart = 0;
  else
    {
      lstart = 1;
      Dev[0] = nullDev;
    }

  /* Path */
  double xwr,xwx,eta,pi,z,v;
  for (int l=lstart;l<L;l++)
    {
      /*Rprintf("l=%d\n",l);*/
      if (l != 0)
	{
	  beta0_old = beta0[l-1];
	  for (int j=0;j<p;j++) beta_old[j] = beta[(l-1)*p+j];
	}

      while (iter[l] < max_iter)
	{
	  converged = 0;
	  iter[l] = iter[l] + 1;

	  /* Check dfmax */
	  active = 0;
	  for (int j=0;j<p;j++) if (beta[l*p+j]!=0) active++;
	  if (active > dfmax)
	    {
	      for (int ll=l;ll<L;ll++)
		{
		  beta0[ll] = R_NaReal;
		  for (int j=0;j<p;j++) beta[ll*p+j] = R_NaReal;
		}
	      free_vector(beta_old);
	      free_vector(w);
	      free_vector(r);
	      return;
	    }

	  /* Approximate L */
	  Dev[l] = 0;
	  for (int i=0;i<n;i++)
	    {
	      eta = beta0_old;
	      for (int j=0;j<p;j++) eta = eta + x[j*n+i]*beta_old[j];
	      pi = exp(eta)/(1+exp(eta));
	      if (pi > .9999)
		{
		  pi = 1;
		  w[i] = .0001;
		}
	      else if (pi < .0001)
		{
		  pi = 0;
		  w[i] = .0001;
		}
	      else w[i] = pi*(1-pi);
	      r[i] = (y[i] - pi)/w[i];
	      Dev[l] = Dev[l] - y[i]*log(pi)-(1-y[i])*log(1-pi);
	      /*yp = yp + pow(y[i]-pi,2);
		yy = yy + pow(y[i]-ybar,2);*/
	    }
	  if (Dev[l]/nullDev < .01)
	    {
	      if (warn) warning("Model saturated; exiting...");
	      for (int ll=l;ll<L;ll++)
		{
		  beta0[ll] = R_NaReal;
		  for (int j=0;j<p;j++) beta[ll*p+j] = R_NaReal;
		}
	      free_vector(beta_old);
	      free_vector(w);
	      free_vector(r);
	      return;
	    }

	  /* Intercept */
	  xwr = xwx = 0;
	  for (int i=0;i<n;i++)
	    {
	      xwr = xwr + w[i]*r[i];
	      xwx = xwx + w[i];
	    }
	  beta0[l] = xwr/xwx + beta0_old;
	  for (int i=0;i<n;i++) r[i] = r[i] - (beta0[l] - beta0_old);

	  /* Covariates */
	  for (int j=0;j<p;j++)
	    {
	      /* Calculate z */
	      xwr=0;
	      xwx=0;
	      for (int i=0;i<n;i++)
		{
		  xwr = xwr + x[j*n+i]*w[i]*r[i];
		  xwx = xwx + x[j*n+i]*w[i]*x[j*n+i];
		}
	      z = xwr/n + (xwx/n)*beta_old[j];
	      v = xwx/n;

	      /* Update beta_j */
	      double l1 = lambda[l] * multiplier[j] * alpha;
	      double l2 = lambda[l] * multiplier[j] * (1-alpha);
	      if (strcmp(penalty,"MCP")==0) beta[l*p+j] = MCP(z, l1, l2, gamma, v);
	      if (strcmp(penalty,"SCAD")==0) beta[l*p+j] = SCAD(z, l1, l2, gamma, v);
	      if (strcmp(penalty,"lasso")==0) beta[l*p+j] = SCAD(z, l1, l2, gamma, v);

	      /* Update r */
	      if (beta[l*p+j] != beta_old[j]) for (int i=0;i<n;i++) r[i] = r[i] - (beta[l*p+j] - beta_old[j])*x[j*n+i];
	    }

	  /* Check for convergence */
	  if (checkConvergence(beta,beta_old,eps,l,p))
	    {
	      converged = 1;
	      break;
	    }
	  beta0_old = beta0[l];
	  for (int j=0;j<p;j++) beta_old[j] = beta[l*p+j];
	}
    }

  free_vector(beta_old);
  free_vector(w);
  free_vector(r);
}

static const R_CMethodDef cMethods[] = {
  {"cdfit_gaussian", (DL_FUNC) &cdfit_gaussian, 17},
  {"cdfit_binomial", (DL_FUNC) &cdfit_binomial, 19},
  NULL
};

void R_init_ncvreg(DllInfo *info)
{
  R_registerRoutines(info,cMethods,NULL,NULL,NULL);
}
