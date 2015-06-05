// [[Rcpp::depends(RcppGSL)]]
#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>

#include "distributions.h"

using namespace Rcpp;

// This contains all distribution-related functions.
// Including, gamma, log-normal, inverse gaussian, positive stable, PVF
// Closed-form Laplace transforms are included when they exist

// The C++ internal functions are of the form (double x, double *p), where p
// is the parameter vector for the distribution
// Exported functions are of the form (NumericVector x, NumericVector p), where
// p is again a parameter vector

NumericVector vectorized_density(NumericVector *x, NumericVector *p, density_fn density) {
  int n = x->size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i] = density((*x)[i], p->begin());
  }
  return out;
}

NumericVector vectorized_deriv_density(NumericVector *x, NumericVector *p, 
                                       deriv_density_fn density, int deriv_idx) {
  int n = x->size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i] = density((*x)[i], p->begin(), deriv_idx);
  }
  return out;
}

///////////////////////////////////////////////////////////////////////////////
// Gamma (Γ)
double dgamma(double x, double *p) {
  double theta = p[0];
  return gsl_ran_gamma_pdf(x, 1/theta, theta);
}

double deriv_dgamma(double x, double *p, int deriv_idx) {
  double theta = p[0];
  double term1, term2, term3;
  term1 = pow((x/theta),(1/theta - 1));
  term2 = exp(-x/theta);
  term3 = log(theta/x) + gsl_sf_psi(1/theta) + x - 1;
  return (term1 * term2 * term3)/(tgamma(1/theta)*pow(theta,3));
}

double lt_dgamma(int p, double s, double* params) {
  double theta = params[0];
  return pow(-1, p) * pow(1/theta, 1/theta) * 
    pow((1/theta) + s, -((1/theta) + p)) * 
    tgamma((1/theta) + p)/tgamma((1/theta));
}

// [[Rcpp::export]]
NumericVector dgamma_c(NumericVector x, NumericVector theta) {
  return vectorized_density(&x, &theta, dgamma);
}

// [[Rcpp::export]]
NumericVector deriv_dgamma_c(NumericVector x, NumericVector theta) {
  return vectorized_deriv_density(&x, &theta, deriv_dgamma, 0);
}

// [[Rcpp::export]]
double lt_dgamma_c(int p, double s, double theta) {
  return lt_dgamma(p, s, &theta);
}
///////////////////////////////////////////////////////////////////////////////
// Log-normal (LN)
double dlognormal(double x, double *p) {
  double theta = p[0];
  return 0;
}

double deriv_dlognormal(double x, double *p, int deriv_idx) {
  double theta = p[0];
  return 0;
}

// [[Rcpp::export]]
NumericVector dlognormal_c(NumericVector x, NumericVector theta) {
  return vectorized_density(&x, &theta, dlognormal);
}

// [[Rcpp::export]]
NumericVector deriv_dlognormal_c(NumericVector x, NumericVector theta) {
  return vectorized_deriv_density(&x, &theta, deriv_dlognormal, 0);
}
///////////////////////////////////////////////////////////////////////////////
// Inverse Gaussian (IG)
double dinvgauss(double x, double *p) {
  double theta = p[0];
  return 0;
}

double deriv_dinvgauss(double x, double *p, int deriv_idx) {
  double theta = p[0];
  return 0;
}

// [[Rcpp::export]]
NumericVector dinvgauss_c(NumericVector x, NumericVector theta) {
  return vectorized_density(&x, &theta, dinvgauss);
}

// [[Rcpp::export]]
NumericVector deriv_dinvgauss_c(NumericVector x, NumericVector theta) {
  return vectorized_deriv_density(&x, &theta, deriv_dinvgauss, 0);
}

///////////////////////////////////////////////////////////////////////////////
// Positive stable (PS)
double dposstab(double x, double *p) {
  double alpha = p[0];
  return 0;
}

double lt_dposstab(int p, double s, double* params) {
  double alpha = params[0];
  return 0;
}

// [[Rcpp::export]]
double lt_dposstab_c(int p, double s, double theta) {
  return lt_dposstab(p, s, &theta);
}
///////////////////////////////////////////////////////////////////////////////
// Power-variance function (PVF)
double dpvf(double x, double *p) {
  double alpha = p[0];
  return 0;
}

double lt_dpvf(int p, double s, double* params) {
  double alpha = params[0];
  return 0;
}

// [[Rcpp::export]]
double lt_dpvf_c(int p, double s, double theta) {
  return lt_dpvf(p, s, &theta);
}
///////////////////////////////////////////////////////////////////////////////
// Log-normal mixture (LNM)

///////////////////////////////////////////////////////////////////////////////
// Pareto-Positive-Stable (PPS)