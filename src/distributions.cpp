// [[Rcpp::depends(RcppGSL)]]
#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_deriv.h>
#include "distributions.h"

using namespace Rcpp;

#define DPOSSTAB_K 100

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

NumericVector vectorized_deriv_deriv_density(NumericVector *x, NumericVector *p, 
                                       deriv_deriv_density_fn density,
                                       int deriv_idx_1, int deriv_idx_2) {
  int n = x->size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i] = density((*x)[i], p->begin(), deriv_idx_1, deriv_idx_2);
  }
  return out;
}

///////////////////////////////////////////////////////////////////////////////
// Numerical derivatives of LT functions
double lt_deriv(double w, void* data) {
  struct lt_params *params = (struct lt_params *)data;
  params->theta[params->deriv_idx] = w;
  return params->lt(params->p, params->s, params->theta);
}

double derivative(double (*f)(double,void*), void* data, double x) {
  gsl_function F;
  double result, abserr;
  
  F.function = f;
  F.params = data;
  
  gsl_deriv_central(&F, x, 1e-8, &result, &abserr);
  
  return result;
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

double deriv_lt_dgamma(int p, double s, double* params, int deriv_idx) {
  // deriv_idx is ignored since there's only one parameter, theta
  double theta = params[0];
  // TODO: maybe clean this up a bit?
  double tmp1 = pow(1/theta, 1/theta);
  double tmp2 = pow(-1, p);
  double tmp3 = tgamma(p + 1/theta);
  double tmp4 = (1/theta + s);
  double tmp5 = (-1/theta - p);
  
  double term1 = tmp1 * tmp2 * tmp3 * pow(tmp4, tmp5) * 
                  (log(tmp4)/pow(theta, 2) - tmp5/(pow(theta, 2) * tmp4) );
    
  double term2 = tmp1 * tmp2 * (-1/pow(theta, 2) - 
                  log(1/theta)/pow(theta, 2)) * tmp3 * pow(tmp4, tmp5);
    
  double term3 = pow(1/theta, 1/theta + 2) * tmp2 * 
                  tmp3 * gsl_sf_psi(p + 1/theta) * pow(tmp4, tmp5);
    
  double term4 = pow(1/theta, 1/theta + 2) * tmp2 * 
                  gsl_sf_psi(1/theta) * tmp3 * pow(tmp4, tmp5);
    
  return (term1 + term2 - term3 + term4)/tgamma(1/theta);
}

double deriv_deriv_lt_dgamma(int p, double s, double* params, int deriv_idx_1, int deriv_idx_2) {
  // deriv_idx is ignored since there's only one parameter, theta
  double theta = params[0];
  
  return
    ((pow((-1),p) * pow((1/theta),(1/theta)) * ((pow(((1/theta) + s),(-((1/theta) + 
      p))) * (log(((1/theta) + s)) * (1/pow(theta,2))) - pow(((1/theta) + 
      s),((-((1/theta) + p)) - 1)) * ((-((1/theta) + p)) * (1/pow(theta,2)))) * 
      (log(((1/theta) + s)) * (1/pow(theta,2))) - pow(((1/theta) + s),(-((1/theta) + 
      p))) * (log(((1/theta) + s)) * (2 * theta/pow((pow(theta,2)),2)) + 1/pow(theta,2)/((1/theta) + 
      s) * (1/pow(theta,2))) - ((pow(((1/theta) + s),((-((1/theta) + p)) - 
      1)) * (log(((1/theta) + s)) * (1/pow(theta,2))) - pow(((1/theta) + 
      s),(((-((1/theta) + p)) - 1) - 1)) * (((-((1/theta) + p)) - 
      1) * (1/pow(theta,2)))) * ((-((1/theta) + p)) * (1/pow(theta,2))) + 
      pow(((1/theta) + s),((-((1/theta) + p)) - 1)) * (1/pow(theta,2) * (1/pow(theta,2)) - 
      (-((1/theta) + p)) * (2 * theta/pow((pow(theta,2)),2))))) - pow((-1),p) * 
      (pow((1/theta),(1/theta)) * (log((1/theta)) * (1/pow(theta,2))) + pow((1/theta),((1/theta) - 
      1)) * ((1/theta) * (1/pow(theta,2)))) * (pow(((1/theta) + s),(-((1/theta) + 
      p))) * (log(((1/theta) + s)) * (1/pow(theta,2))) - pow(((1/theta) + 
      s),((-((1/theta) + p)) - 1)) * ((-((1/theta) + p)) * (1/pow(theta,2)))) - 
      (pow((-1),p) * (pow((1/theta),(1/theta)) * (log((1/theta)) * (1/pow(theta,2))) + 
      pow((1/theta),((1/theta) - 1)) * ((1/theta) * (1/pow(theta,2)))) * 
      (pow(((1/theta) + s),(-((1/theta) + p))) * (log(((1/theta) + 
      s)) * (1/pow(theta,2))) - pow(((1/theta) + s),((-((1/theta) + 
      p)) - 1)) * ((-((1/theta) + p)) * (1/pow(theta,2)))) - 
      pow((-1),p) * (pow((1/theta),((1/theta) - 1)) * ((1/theta) * (2 * 
      theta/pow((pow(theta,2)),2)) + 1/pow(theta,2) * (1/pow(theta,2))) + (pow((1/theta),((1/theta) - 
      1)) * (log((1/theta)) * (1/pow(theta,2))) + pow((1/theta),(((1/theta) - 
      1) - 1)) * (((1/theta) - 1) * (1/pow(theta,2)))) * ((1/theta) * 
      (1/pow(theta,2))) + (pow((1/theta),(1/theta)) * (log((1/theta)) * 
      (2 * theta/pow((pow(theta,2)),2)) + 1/pow(theta,2)/(1/theta) * (1/pow(theta,2))) + 
      (pow((1/theta),(1/theta)) * (log((1/theta)) * (1/pow(theta,2))) + 
      pow((1/theta),((1/theta) - 1)) * ((1/theta) * (1/pow(theta,2)))) * 
      (log((1/theta)) * (1/pow(theta,2))))) * pow(((1/theta) + 
      s),(-((1/theta) + p))))) * tgamma((1/theta) + p) - 
      (pow((-1),p) * pow((1/theta),(1/theta)) * (pow(((1/theta) + s),(-((1/theta) + 
      p))) * (log(((1/theta) + s)) * (1/pow(theta,2))) - pow(((1/theta) + 
      s),((-((1/theta) + p)) - 1)) * ((-((1/theta) + p)) * (1/pow(theta,2)))) - 
      pow((-1),p) * (pow((1/theta),(1/theta)) * (log((1/theta)) * (1/pow(theta,2))) + 
      pow((1/theta),((1/theta) - 1)) * ((1/theta) * (1/pow(theta,2)))) * 
      pow(((1/theta) + s),(-((1/theta) + p)))) * (1/pow(theta,2) * 
      (tgamma((1/theta) + p) * gsl_sf_psi((1/theta) + p))) - ((pow((-1),p) * 
      pow((1/theta),(1/theta)) * (pow(((1/theta) + s),(-((1/theta) + p))) * 
      (log(((1/theta) + s)) * (1/pow(theta,2))) - pow(((1/theta) + s),((-((1/theta) + 
      p)) - 1)) * ((-((1/theta) + p)) * (1/pow(theta,2)))) - pow((-1),p) * 
      (pow((1/theta),(1/theta)) * (log((1/theta)) * (1/pow(theta,2))) + pow((1/theta),((1/theta) - 
      1)) * ((1/theta) * (1/pow(theta,2)))) * pow(((1/theta) + s),(-((1/theta) + 
      p)))) * (1/pow(theta,2) * (tgamma((1/theta) + p) * gsl_sf_psi((1/theta) + 
      p))) - pow((-1),p) * pow((1/theta),(1/theta)) * pow(((1/theta) + s),(-((1/theta) + 
      p))) * (1/pow(theta,2) * (tgamma((1/theta) + p) * (1/pow(theta,2) * gsl_sf_psi_1((1/theta) + 
      p)) + 1/pow(theta,2) * (tgamma((1/theta) + p) * gsl_sf_psi((1/theta) + 
      p)) * gsl_sf_psi((1/theta) + p)) + 2 * theta/pow((pow(theta,2)),2) * (tgamma((1/theta) + 
      p) * gsl_sf_psi((1/theta) + p)))))/tgamma((1/theta)) + ((pow((-1),p) * 
      pow((1/theta),(1/theta)) * (pow(((1/theta) + s),(-((1/theta) + p))) * 
      (log(((1/theta) + s)) * (1/pow(theta,2))) - pow(((1/theta) + s),((-((1/theta) + 
      p)) - 1)) * ((-((1/theta) + p)) * (1/pow(theta,2)))) - pow((-1),p) * 
      (pow((1/theta),(1/theta)) * (log((1/theta)) * (1/pow(theta,2))) + pow((1/theta),((1/theta) - 
      1)) * ((1/theta) * (1/pow(theta,2)))) * pow(((1/theta) + s),(-((1/theta) + 
      p)))) * tgamma((1/theta) + p) - pow((-1),p) * pow((1/theta),(1/theta)) * 
      pow(((1/theta) + s),(-((1/theta) + p))) * (1/pow(theta,2) * (tgamma((1/theta) + 
      p) * gsl_sf_psi((1/theta) + p)))) * (1/pow(theta,2) * (tgamma((1/theta)) * 
      gsl_sf_psi((1/theta))))/pow(tgamma((1/theta)),2) + ((((pow((-1),p) * pow((1/theta),(1/theta)) * 
      (pow(((1/theta) + s),(-((1/theta) + p))) * (log(((1/theta) + s)) * 
      (1/pow(theta,2))) - pow(((1/theta) + s),((-((1/theta) + p)) - 
      1)) * ((-((1/theta) + p)) * (1/pow(theta,2)))) - pow((-1),p) * (pow((1/theta),(1/theta)) * 
      (log((1/theta)) * (1/pow(theta,2))) + pow((1/theta),((1/theta) - 1)) * 
      ((1/theta) * (1/pow(theta,2)))) * pow(((1/theta) + s),(-((1/theta) + 
      p)))) * tgamma((1/theta) + p) - pow((-1),p) * pow((1/theta),(1/theta)) * 
      pow(((1/theta) + s),(-((1/theta) + p))) * (1/pow(theta,2) * (tgamma((1/theta) + 
      p) * gsl_sf_psi((1/theta) + p)))) * (1/pow(theta,2) * (tgamma((1/theta)) * 
      gsl_sf_psi((1/theta)))) - pow((-1),p) * pow((1/theta),(1/theta)) * pow(((1/theta) + 
      s),(-((1/theta) + p))) * tgamma((1/theta) + p) * (1/pow(theta,2) * 
      (tgamma((1/theta)) * (1/pow(theta,2) * gsl_sf_psi_1((1/theta))) + 1/pow(theta,2) * 
      (tgamma((1/theta)) * gsl_sf_psi((1/theta))) * gsl_sf_psi((1/theta))) + 
      2 * theta/pow((pow(theta,2)),2) * (tgamma((1/theta)) * gsl_sf_psi((1/theta)))))/pow(tgamma((1/theta)),2) + 
      pow((-1),p) * pow((1/theta),(1/theta)) * pow(((1/theta) + s),(-((1/theta) + 
      p))) * tgamma((1/theta) + p) * (1/pow(theta,2) * (tgamma((1/theta)) * 
      gsl_sf_psi((1/theta)))) * (2 * (1/pow(theta,2) * (tgamma((1/theta)) * 
      gsl_sf_psi((1/theta))) * tgamma((1/theta))))/pow((pow(tgamma((1/theta)),2)),2));
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

// [[Rcpp::export]]
double deriv_lt_dgamma_c(int p, double s, double theta) {
  return deriv_lt_dgamma(p, s, &theta, 0);
}

// [[Rcpp::export]]
double deriv_deriv_lt_dgamma_c(int p, double s, double theta) {
  return deriv_deriv_lt_dgamma(p, s, &theta, 0, 0);
}

///////////////////////////////////////////////////////////////////////////////
// Log-normal (LN)
double dlognormal(double x, double *params) {
  double theta = params[0];
  double factor1 = 1/(x*sqrt(theta*2*PI));
  double factor2 = exp(-pow(log(x),2)/(2*theta));
  return factor1*factor2;
}

double deriv_dlognormal(double x, double *params, int deriv_idx) {
  double theta = params[0];
  double term1_numer = pow(log(x),2) * exp(-pow(log(x),2)/(2*theta));
  double term1_denom = 2 * sqrt(2*PI) * pow(theta, 5.0/2.0) * x;
  double term2_numer = exp(-pow(log(x),2)/(2*theta));
  double term2_denom = 2 * sqrt(2*PI) * pow(theta, 3.0/2.0) * x;
  return (term1_numer/term1_denom) - (term2_numer/term2_denom);
}

double deriv_deriv_dlognormal(double x, double *params, int deriv_idx_1, int deriv_idx_2) {
  double theta = params[0];
  return pow(log(x), 2) * (exp(-(pow(log(x), 2))/(2 * theta)) * ((pow(log(x), 2)) * 
   2/pow((2 * theta),2)))/(2 * sqrt(2 * PI) * pow(theta,5.0/2.0) * x) - (pow(log(x), 2) * 
   exp(-(pow(log(x), 2))/(2 * theta))) * (2 * sqrt(2 * PI) * 
   (pow(theta,((5.0/2.0) - 1)) * (5/2)) * x)/pow((2 * sqrt(2 * PI) * pow(theta,5.0/2.0) * x),2) - 
   (exp(-(pow(log(x), 2))/(2 * theta)) * ((pow(log(x), 2)) * 2/pow((2 * theta),2))/(2 * 
   sqrt(2 * PI) * pow(theta,3.0/2.0) * x) - (exp(-(pow(log(x), 2))/(2 * 
   theta))) * (2 * sqrt(2 * PI) * (pow(theta,((3.0/2.0) - 1)) * (3.0/2.0)) * 
   x)/pow((2 * sqrt(2 * PI) * pow(theta,3.0/2.0) * x),2));
}

// [[Rcpp::export]]
NumericVector dlognormal_c(NumericVector x, NumericVector theta) {
  return vectorized_density(&x, &theta, dlognormal);
}

// [[Rcpp::export]]
NumericVector deriv_dlognormal_c(NumericVector x, NumericVector theta) {
  return vectorized_deriv_density(&x, &theta, deriv_dlognormal, 0);
}

// [[Rcpp::export]]
NumericVector deriv_deriv_dlognormal_c(NumericVector x, NumericVector theta) {
  return vectorized_deriv_deriv_density(&x, &theta, deriv_deriv_dlognormal, 0, 0);
}

///////////////////////////////////////////////////////////////////////////////
// Inverse Gaussian (IG)
double dinvgauss(double x, double *p) {
  double theta = p[0];
  double numer = exp(-pow(x - 1, 2)/(2*theta*x));
  double denom = sqrt(2*PI*theta*pow(x,3));
  return numer/denom;
}

double deriv_dinvgauss(double x, double *p, int deriv_idx) {
  double theta = p[0];
  double term1_numer = pow(x - 1, 2) * exp(-pow(x - 1, 2)/(2*theta*x));
  double term1_denom = 2 * sqrt(2*PI) * pow(theta, 2) * x * sqrt(theta * pow(x, 3));
  double term2_numer = pow(x, 3) * exp(-pow(x - 1, 2)/(2*theta*x));
  double term2_denom = 2 * sqrt(2*PI) * pow(theta * pow(x, 3), 3.0/2.0);
  return term1_numer/term1_denom - term2_numer/term2_denom;
}

double deriv_deriv_dinvgauss(double x, double *p, int deriv_idx_1, int deriv_idx_2) {
  double theta = p[0];
  return (exp(-(pow((x - 1),2))/(2 * theta * x)) * ((pow((x - 1),2)) * (2 * x)/pow(2 * 
          theta * x,2)) * ((pow((x - 1),2)) * (2 * x)/pow(2 * theta * x,2)) - 
          exp(-(pow((x - 1),2))/(2 * theta * x)) * ((pow((x - 1),2)) * (2 * x) * 
          (2 * (2 * x * (2 * theta * x)))/pow(pow(2 * theta * x,2),2)))/sqrt(2 * 
          PI * theta * (pow(x,3))) - exp(-(pow((x - 1),2))/(2 * theta * x)) * 
          ((pow((x - 1),2)) * (2 * x)/pow(2 * theta * x,2)) * (0.5 * (2 * PI * 
          (pow(x,3)) * pow(2 * PI * theta * (pow(x,3)),-0.5)))/pow(sqrt(2 * PI * theta * 
          (pow(x,3))),2) - ((exp(-(pow((x - 1),2))/(2 * theta * x)) * ((pow((x - 1),2)) * 
          (2 * x)/pow(2 * theta * x,2)) * (0.5 * (2 * PI * (pow(x,3)) * pow(2 * 
          PI * theta * (pow(x,3)),-0.5))) - exp(-(pow((x - 1),2))/(2 * theta * 
          x)) * (0.5 * (2 * PI * (pow(x,3)) * (pow(2 * PI * theta * (pow(x,3)),-(0.5 + 
          1)) * (0.5 * (2 * PI * (pow(x,3))))))))/pow(sqrt(2 * PI * theta * (pow(x,3))),2) - 
          exp(-(pow((x - 1),2))/(2 * theta * x)) * (0.5 * (2 * PI * (pow(x,3)) * 
          pow(2 * PI * theta * (pow(x,3)),-0.5))) * (2 * (0.5 * (2 * PI * 
          (pow(x,3)) * pow((2 * PI * theta * (pow(x,3))),-0.5)) * sqrt(2 * PI * 
          theta * (pow(x,3)))))/pow(pow(sqrt(2 * PI * theta * (pow(x,3))),2),2));
}

// [[Rcpp::export]]
NumericVector dinvgauss_c(NumericVector x, NumericVector theta) {
  return vectorized_density(&x, &theta, dinvgauss);
}

// [[Rcpp::export]]
NumericVector deriv_dinvgauss_c(NumericVector x, NumericVector theta) {
  return vectorized_deriv_density(&x, &theta, deriv_dinvgauss, 0);
}

// [[Rcpp::export]]
NumericVector deriv_deriv_dinvgauss_c(NumericVector x, NumericVector theta) {
  return vectorized_deriv_deriv_density(&x, &theta, deriv_deriv_dinvgauss, 0, 0);
}

///////////////////////////////////////////////////////////////////////////////
// Power-variance function (PVF)
double dposstab(double x, double *params) {
  double alpha = params[0];
  
  double factor1 = -1/(PI * x);
  double factor2 = 0;
  for (int k = 1; k <= DPOSSTAB_K; ++k) {
    factor2 += tgamma(k * alpha + 1) / gsl_sf_fact(k) * pow(-pow(x, -alpha), k) * sin(alpha * k * PI);
  }
  
  return factor1*factor2;
}

double dpvf(double x, double *params) {
  double alpha = params[0];
  double factor1 = dposstab(x * pow(alpha, 1/alpha), params);
  double factor2 = pow(alpha, 1/alpha) * exp(-x) * exp(1/alpha);
  return factor1*factor2;
}

// LT coefficients for PS and PVF distributions
// [[Rcpp::export]]
double lt_dpvf_coef(int p, int j, double alpha) {
  // Base cases
  if (p == j) return 1;
  if (j == 1) return (tgamma(p - alpha)/tgamma(1 - alpha));
      
  double term1 = lt_dpvf_coef(p - 1, j - 1, alpha);
  double term2 = lt_dpvf_coef(p - 1, j, alpha)*((p - 1) - j*alpha);
  return term1 + term2;
}

// [[Rcpp::export]]
double lt_deriv_dpvf_coef(int p, int j, double alpha) {
  // Base cases
  if (p == j) return 0;
  if (j == 1) {
    double numer = tgamma(p - alpha)*(gsl_sf_psi(1 - alpha) - gsl_sf_psi(p - alpha));
    double denom = tgamma(1 - alpha);
    return numer/denom;
  }
    
  double term1 = lt_deriv_dpvf_coef(p - 1, j - 1, alpha);
  double term2 = lt_deriv_dpvf_coef(p - 1, j, alpha) * ((p - 1) - j*alpha);
  double term3 = -j * lt_dpvf_coef(p - 1, j, alpha);
  return term1 + term2 + term3;
}

// [[Rcpp::export]]
double lt_deriv_deriv_dpvf_coef(int m, int j, double alpha) {
  // Base cases
  if (m == j) return 0;
    
  if (j == 1) {
    return((tgamma(m - alpha) * gsl_sf_psi_1(m - alpha) + tgamma(m - alpha) * 
           gsl_sf_psi(m - alpha) * gsl_sf_psi(m - alpha))/tgamma(1 - alpha) - 
           tgamma(m - alpha) * gsl_sf_psi(m - alpha) * (tgamma(1 - alpha) * 
           gsl_sf_psi(1 - alpha))/pow(tgamma(1 - alpha),2) - ((tgamma(m - 
           alpha) * (tgamma(1 - alpha) * gsl_sf_psi_1(1 - alpha) + tgamma(1 - 
           alpha) * gsl_sf_psi(1 - alpha) * gsl_sf_psi(1 - alpha)) + tgamma(m - 
           alpha) * gsl_sf_psi(m - alpha) * (tgamma(1 - alpha) * gsl_sf_psi(1 - 
           alpha)))/pow(tgamma(1 - alpha),2) - tgamma(m - alpha) * (tgamma(1 - 
           alpha) * gsl_sf_psi(1 - alpha)) * (2 * (tgamma(1 - alpha) * gsl_sf_psi(1 - 
           alpha) * tgamma(1 - alpha)))/pow((pow(tgamma(1 - alpha),2)),2)));
  }
    
  double term1 = lt_deriv_deriv_dpvf_coef(m - 1, j - 1, alpha);
  double term2 = lt_deriv_deriv_dpvf_coef(m - 1, j, alpha) * ((m - 1) - j*alpha);
  double term3 = -j * lt_deriv_dpvf_coef(m - 1, j, alpha);
      
  return term1 + term2 + 2*term3;
}

double lt_dpvf(int p, double s, double* params) {
  double alpha = params[0];
  
  if (p == 0) {
    return exp(-(pow(1 + s, alpha) - 1)/alpha);
  }
  
  double factor1 = lt_dpvf(0, s, params);
  
  double factor2 = 0;
  for (int j = 1; j <= p; ++j) {
    factor2 += lt_dpvf_coef(p, j, alpha) * pow(1 + s, j*alpha - p);
  }
  
  return pow(-1, p)*factor1*factor2;
}

double deriv_lt_dpvf(int p, double s, double* params, int deriv_idx) {
  // deriv_idx ignored since only one parameter
  double alpha = params[0];
  
  if (p == 0) {
    double factor1 = exp((1 - pow(s + 1, alpha))/alpha);
    double factor2_term1 = -(1 - pow(s + 1, alpha))/pow(alpha, 2);
    double factor2_term2 = -(pow(s + 1, alpha) * log(s + 1))/alpha;
    return factor1*(factor2_term1 + factor2_term2);
  }
  
  lt_params lt_p = (struct lt_params){p, s, params, lt_dpvf, deriv_idx};
  return derivative(&lt_deriv, &lt_p, params[deriv_idx]);
}

double deriv_deriv_lt_dpvf(int m, double s, double* params, int deriv_idx_1, int deriv_idx_2) {
  // deriv_idx ignored since only one parameter
  double alpha = params[0];
  
  if (m == 0) {
    return -(exp(-(pow((1 + s),alpha) - 1)/alpha) * (pow((1 + s),alpha) * log((1 + 
             s)) * log((1 + s))/alpha - pow((1 + s),alpha) * log((1 + s))/pow(alpha,2) - 
             (pow((1 + s),alpha) * log((1 + s))/pow(alpha,2) - (pow((1 + s),alpha) - 
             1) * (2 * alpha)/pow((pow(alpha,2)),2))) - exp(-(pow((1 + s),alpha) - 
             1)/alpha) * (pow((1 + s),alpha) * log((1 + s))/alpha - (pow((1 + s),alpha) - 
             1)/pow(alpha,2)) * (pow((1 + s),alpha) * log((1 + s))/alpha - (pow((1 + 
             s),alpha) - 1)/pow(alpha,2)));
  }
  
  double sum1, sum2, sum3, sum4, sum5, sum6;
  sum1 = sum2 = sum3 = sum4 = sum5 = sum6 = 0;
  
  for (int j = 1; j <= m; ++j) {
    sum1 += lt_dpvf_coef(m, j, alpha) * pow((1 + s),(j*alpha - m)) * pow(j,2) * pow(log(1 + s),2);
    sum2 += lt_deriv_dpvf_coef(m, j, alpha) * pow((1 + s),(j*alpha - m)) * j * log(1 + s);
    sum3 += lt_deriv_deriv_dpvf_coef(m, j, alpha) * pow((1 + s),(j*alpha - m));
    sum4 += lt_dpvf_coef(m, j, alpha) * pow((1 + s),(j*alpha - m)) * j * log(1 + s);
    sum5 += lt_dpvf_coef(m, j, alpha) * pow((1 + s),(j*alpha - m));
    sum6 += lt_deriv_dpvf_coef(m, j, alpha) * pow((1 + s),(j*alpha - m));
  }
  
  double lt = lt_dpvf(0, s, params);
  double dlt = deriv_lt_dpvf(0, s, params, deriv_idx_1);
  double ddlt = deriv_deriv_lt_dpvf(0, s, params, deriv_idx_1, deriv_idx_2);
  
  double term1 = pow(-1, m) * lt * sum1;
  double term2 = 2 * pow(-1, m) * lt * sum2;
  double term3 = pow(-1, m) * lt * sum3;
  double term4 = 2 * pow(-1, m) * dlt * sum4;
  double term5 = pow(-1, m) * ddlt * sum5;
  double term6 = 2 * pow(-1, m) * dlt * sum6;
    
  return term1 + term2 + term3 + term4 + term5 + term6;
}

// [[Rcpp::export]]
NumericVector dpvf_c(NumericVector x, NumericVector alpha) {
  return vectorized_density(&x, &alpha, dpvf);
}

// [[Rcpp::export]]
double lt_dpvf_c(int m, double s, double alpha) {
  return lt_dpvf(m, s, &alpha);
}

// [[Rcpp::export]]
double deriv_lt_dpvf_c(int m, double s, double alpha) {
  return deriv_lt_dpvf(m, s, &alpha, 0);
}

// [[Rcpp::export]]
double deriv_deriv_lt_dpvf_c(int m, double s, double alpha) {
  return deriv_deriv_lt_dpvf(m, s, &alpha, 0, 0);
}

///////////////////////////////////////////////////////////////////////////////
// Positive stable (PS)

// Laplace transform of positive stable
double lt_dposstab(int p, double s, double* params) {
  double alpha = params[0];
  
  if (p == 0) {
    return exp(-alpha*pow(s, alpha)/alpha);
  }
  
  double factor1 = lt_dposstab(0, s, params);
  
  double factor2 = 0;
  for (int j = 1; j <= p; ++j) {
    factor2 += lt_dpvf_coef(p, j, alpha) * pow(alpha, j) * pow(s, (j*alpha - p));
  }
  
  return pow(-1, p)*factor1*factor2;
}

// Derivative wrt. alpha of the Laplace transform of positive stable
double deriv_lt_dposstab(int p, double s, double* params, int deriv_idx) {
  // deriv_idx ignored since only one parameter
  double alpha = params[0];
  
  if (p == 0) {
    return -exp(-pow(s, alpha)) * pow(s, alpha) * log(s);
  }
  
  lt_params lt_p = (struct lt_params){p, s, params, lt_dposstab, deriv_idx};
  return derivative(&lt_deriv, &lt_p, params[deriv_idx]);
}

// [[Rcpp::export]]
NumericVector dposstab_c(NumericVector x, NumericVector alpha) {
  return vectorized_density(&x, &alpha, dposstab);
}

// [[Rcpp::export]]
double lt_dposstab_c(int p, double s, double alpha) {
  return lt_dposstab(p, s, &alpha);
}

// [[Rcpp::export]]
double deriv_lt_dposstab_c(int p, double s, double alpha) {
  return deriv_lt_dposstab(p, s, &alpha, 0);
}

///////////////////////////////////////////////////////////////////////////////
// Log-normal mixture (LNM)

///////////////////////////////////////////////////////////////////////////////
// Pareto-Positive-Stable (PPS)

///////////////////////////////////////////////////////////////////////////////
// Misc. functions

// Reimann zeta function, modified from package degreenet
// [[Rcpp::export]]
double zeta(double s) {
  int a, k;
  int m, n, m2;
  double p, a2, fred, sum;
  
  double b2[12];
  b2[0] = 1.0 / 6.0;
  b2[1] = -1.0 / 30.0;
  b2[2] = 1.0 / 42.0;
  b2[3] = -1.0 / 30.0;
  b2[4] = 5.0 / 66.0;
  b2[5] = -691.0 / 2730.0;
  b2[6] = 7.0 / 6.0;
  b2[7] = -3617.0 / 510.0;
  b2[8] = 4386.7 / 79.8;
  b2[9] = -1746.11 / 3.30;
  b2[10] = 8545.13 / 1.38;
  b2[11] = -2363.64091 / 0.02730;
  a = 12;
  k = 8;
  a2 = a * a;
  p = s / 2.0 / a2;
  sum = 1.0 / (s - 1.0) + 0.5 / a + b2[0] * p;
  for (m = 1; m < k; m++){
    m2 = m + m + 2;
    p=p*(s+m2-3.0)*(s+m2-2.0)/(m2-1.0)/m2/a2;
    sum = sum + p * b2[m];
  }
  fred = exp((s - 1.0) * log(1.0 * a));
  sum = 1.0 + sum / fred;
  for (n = 2; n < a; n++){
    sum = sum + 1.0 * exp(-s * log(1.0 * n));
  }
  return sum;
}