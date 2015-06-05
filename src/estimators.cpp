// [[Rcpp::depends(RcppGSL)]]
#include <cmath>
#include <RcppGSL.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>
#include "distributions.h"
using namespace Rcpp;

// Args to the adaptive quadrature integration
// integrate from (LOWER, infty), which is mapped to (0, 1] by  x = a + (1-t)/t
#define LOWER 0 

// Convergence of integration reached when either the absolute or relative error 
// is within these tolerances
#define EPSABS 0
#define EPSREL 1e-8

// Size of the work space and max number of subintervals
#define LIMIT 1e8

// The integration wrapper
double integrate(double (*f)(double,void*), void* data) {
  gsl_integration_workspace *work = gsl_integration_workspace_alloc(LIMIT);
  double result, error;
  
  gsl_function F;
  F.function = f;
  F.params = data;
  
  int ret = gsl_integration_qagiu(&F, LOWER, EPSABS, EPSREL, LIMIT, work, &result, &error);
  gsl_integration_workspace_free(work);
  return result;
}

// derivative of phi_ki(gamma, lambda, t), used for numerical integration
double phi_deriv(double w, void* data) {
  struct phi_params *params = (struct phi_params *)data;
  return pow(w, params->N_dot + params->k - 1) * exp(-w * params->H_dot) * 
    params->density(w, params->theta);
}

// derivative of phi_ki(gamma, lambda, t), where f' is used
double phi_prime_deriv(double w, void* data) {
  struct phi_prime_params *params = (struct phi_prime_params *)data;
  return pow(w, params->N_dot + params->k - 1) * exp(-w * params->H_dot) * 
    params->deriv_density(w, params->theta, params->deriv_idx);
}

// The following phi functions are the only places that branch based on the
// frailty distribution. If a new distribution is implemented, it needs to go
// here, with the respective density and deriv density functions (or LT)

double phi(int k, int N_dot, double H_dot, double *theta, String frailty) {
  // Laplace transform integrals
  if (frailty == "gamma") {
    return lt_dgamma(N_dot + k - 1, H_dot, theta) * pow(-1, N_dot + k - 1);
    // phi_p = (struct phi_params){k, N_dot, H_dot, theta, dgamma};
  } else if (frailty == "posstab") {
    return lt_dposstab(N_dot + k - 1, H_dot, theta) * pow(-1, N_dot + k - 1);
  } else if (frailty == "pvf") {
    return lt_dpvf(N_dot + k - 1, H_dot, theta) * pow(-1, N_dot + k - 1);
  }
  
  // Numerical integration
  phi_params phi_p;
  if (frailty == "lognormal") {
    phi_p = (struct phi_params){k, N_dot, H_dot, theta, dlognormal};
  } if (frailty == "invgauss") {
    phi_p = (struct phi_params){k, N_dot, H_dot, theta, dinvgauss};
  } else {
    throw std::range_error("Unsupported frailty distribution");
  }
  
  return integrate(&phi_deriv, &phi_p);
}

// phi using the derivative of the density wrt. parameter p[derive_idx]
double phi_prime(int k, int N_dot, double H_dot, double *theta, String frailty, int deriv_idx) {
  // TODO: Laplace transform integrals?
  if (frailty == "posstab") {
    throw std::range_error("Unsupported frailty distribution");
  } else if (frailty == "pvf") {
    throw std::range_error("Unsupported frailty distribution");
  }
  
  // Numerical integration
  phi_prime_params phi_p;
  if (frailty == "gamma") {
    phi_p = (struct phi_prime_params){k, N_dot, H_dot, theta, deriv_dgamma, deriv_idx};
  } else if (frailty == "lognormal") {
    phi_p = (struct phi_prime_params){k, N_dot, H_dot, theta, deriv_dlognormal, deriv_idx};
  } else if (frailty == "invgauss") {
    phi_p = (struct phi_prime_params){k, N_dot, H_dot, theta, deriv_dinvgauss, deriv_idx};
  } else {
    throw std::range_error("Unsupported frailty distribution");
  }
  
  return integrate(&phi_prime_deriv, &phi_p);
}

double psi(int N_dot, double H_dot, double* theta, String frailty) {
  // Shortcut for gamma frailty
//   if (frailty == "gamma") {
//     return (N_dot + 1/theta[0])/(H_dot + 1/theta[0]);
//   }
  return phi(2, N_dot, H_dot, theta, frailty)/phi(1, N_dot, H_dot, theta, frailty);
}

// [[Rcpp::export]]
Rcpp::List baseline_hazard_estimator(Rcpp::List X_, Rcpp::List k_,
                                     NumericVector d_,
                                     Rcpp::List Y_, Rcpp::List N_dot,
                                     NumericVector beta, NumericVector theta, 
                                     String frailty) {
  int n_timesteps = d_.size();
  int n_clusters = X_.size();
  
  // H_ is a list of matrices
  Rcpp::List H_ = clone(Y_); //Rcpp::List::create();
  
  // H_dot is a list of numeric vectors
  Rcpp::List H_dot = clone(N_dot);
  
  for (int i = 0; i < n_clusters; ++i) {
    // H_[i](NumericMatrix(n_members, n_timesteps));
    Rcpp::NumericMatrix H_i = H_[i];
    Rcpp::NumericVector H_dot_i = H_dot[i];
    for (int k = 0; k < n_timesteps; ++k) {
      for (int j = 0; j < H_i.nrow(); ++j) {
        H_i(j, k) = 0;
      }
      H_dot_i[k] = 0;
    }
  }
  // lambda_hat is the baseline cumulative hazard estimate
  NumericVector lambda_hat(n_timesteps);
  NumericVector delta_lambda_hat(n_timesteps);
  
  double denom;
  // k starts from 1, not a typo
  for (int k = 1; k < n_timesteps; ++k) {
    denom = 0;
    for (int i = 0; i < n_clusters; ++i) {
      Rcpp::NumericMatrix X_i = X_[i];
      
      Rcpp::NumericMatrix Y_i = Y_[i];
      Rcpp::NumericVector N_dot_i = N_dot[i];
      Rcpp::NumericVector H_dot_i = H_dot[i];
      
      double tmp = 0;
      for (int j = 0; j < X_i.nrow(); ++j) {
        tmp += Y_i(j, k) * exp(sum(beta * X_i(j, _)));
      }
      
      denom += tmp * psi(N_dot_i[k-1], H_dot_i[k-1], theta.begin(), frailty);
    }
    
    delta_lambda_hat[k] = d_[k]/denom;
    
    lambda_hat[k] = lambda_hat[k-1] + delta_lambda_hat[k];
    
    for (int i = 0; i < n_clusters; ++i) {
      Rcpp::NumericMatrix X_i = X_[i];
      Rcpp::NumericVector k_i = k_[i];
      
      Rcpp::NumericMatrix H_i = H_[i];
      Rcpp::NumericVector H_dot_i = H_dot[i];
      
      for (int j = 0; j < X_i.nrow(); ++j) {
        // k_i[j] - 1 since k_ is the rank of failure as an R index
        int k_min = min(NumericVector::create(k_i[j] - 1, k));
        H_i(j, k) = lambda_hat[k_min] * exp(sum(beta * X_i(j, _)));
      }
      
      H_dot_i[k] = sum(H_i( _ , k));
    }
  }
  return Rcpp::List::create(Rcpp::Named("H_") = H_,
                            Rcpp::Named("H_dot") = H_dot,
                            Rcpp::Named("lambda_hat") = lambda_hat);
}

// [[Rcpp::export]]
double U_r(Rcpp::List X_,
           Rcpp::List R_, 
           Rcpp::List I_, 
           Rcpp::List N_dot,
           Rcpp::List H_, 
           Rcpp::List H_dot,
           NumericVector beta, 
           NumericVector theta,
           int beta_idx, String frailty) {
  
  int n_clusters = X_.size();
  
  double term1 = 0;
  double term2 = 0;
  for (int i = 0; i < n_clusters; ++i) {
    Rcpp::NumericMatrix X_i = X_(i);
    Rcpp::NumericMatrix H_i = H_(i);
    Rcpp::NumericVector I_i = I_(i); // Failure indicator
    Rcpp::NumericVector R_i = R_(i); // Failure rank, whare T[r] = failure time
    Rcpp::NumericVector N_dot_i = N_dot(i);
    Rcpp::NumericVector H_dot_i = H_dot(i);
    int tau_k = N_dot_i.size() - 1;
    
    double tmp = 0;
    for (int j = 0; j < X_i.nrow(); j++) {
      term1 += I_i[j] * X_i(j,beta_idx - 1);
      tmp += H_i(j, R_i[j] - 1) * X_i(j,beta_idx - 1);
    }
    term2 += tmp * psi(N_dot_i(tau_k), 
                       H_dot_i(tau_k), 
                       theta.begin(), "gamma");
  }
  term1 = term1/n_clusters;
  term2 = term2/n_clusters;
  return term1 - term2;
}

// [[Rcpp::export]]
double U_p(Rcpp::List X_,
           Rcpp::List R_, 
           Rcpp::List I_, 
           Rcpp::List N_dot,
           Rcpp::List H_, 
           Rcpp::List H_dot,
           NumericVector beta, 
           NumericVector theta,
           int theta_idx, String frailty) {
  
  int n_clusters = X_.size();
  
  double out = 0;
  for (int i = 0; i < n_clusters; ++i) {
    Rcpp::NumericMatrix X_i = X_(i);
    Rcpp::NumericMatrix H_i = H_(i);
    Rcpp::NumericVector I_i = I_(i); // Failure indicator
    Rcpp::NumericVector R_i = R_(i);
    Rcpp::NumericVector N_dot_i = N_dot(i);
    Rcpp::NumericVector H_dot_i = H_dot(i);
    int tau_k = N_dot_i.size() - 1;
    
    double denom = phi(1, 
                       N_dot_i(tau_k), 
                       H_dot_i(tau_k),
                       theta.begin(), 
                       frailty);
    
    double numer = phi_prime(1, 
                             N_dot_i(tau_k),
                             H_dot_i(tau_k),
                             theta.begin(), 
                             frailty, 
                             theta_idx - 1); // theta_idx is an R index
    out += numer/denom;
  }
  
  out = out/n_clusters;
  
  return out;
}

// [[Rcpp::export]]
double log_likelihood(List X_, 
                      List R_, 
                      List I_, 
                      List N_dot,
                      List H_dot,
                      NumericVector lambda, 
                      NumericVector beta, 
                      NumericVector theta, 
                      String frailty) {
  
  int n_clusters = X_.size();
  
  double term1 = 0;
  double term2 = 0;
  for (int i = 0; i < n_clusters; ++i) {
    Rcpp::NumericMatrix X_i = X_(i);
    Rcpp::NumericVector I_i = I_(i); // Failure indicator
    Rcpp::NumericVector R_i = R_(i); // Failure rank, whare T[r] = failure time
    Rcpp::NumericVector N_dot_i = N_dot(i);
    Rcpp::NumericVector H_dot_i = H_dot(i);
    int tau_k = N_dot_i.size() - 1;
    
    double term1_tmp = 0;
    for (int j = 0; j < X_i.nrow(); j++) {
      // R_i[i] is an R index
      term1_tmp += I_i(j) * log(lambda(R_i(j) - 1) * exp(sum(beta * X_i(j, _))));
    }
    
    term1 += term1_tmp;
    
    term2 += log(phi(1, 
                     N_dot_i(tau_k), 
                     H_dot_i(tau_k), 
                     theta.begin(), frailty));
  }
  
  return term1 + term2;
}