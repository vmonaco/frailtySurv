// [[Rcpp::depends(RcppGSL)]]
#include <cmath>
#include <RcppGSL.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>
#include "distributions.h"

using namespace Rcpp;

// Args to the adaptive quadrature integration
// integrate from (INT_LOWER, infty), which is mapped to (0, 1] by  x = a + (1-t)/t
#define INT_LOWER 0 

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
  
  int ret = gsl_integration_qagiu(&F, INT_LOWER, EPSABS, EPSREL, LIMIT, work, &result, &error);
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

// derivative of phi_ki(gamma, lambda, t), where f' is used
double phi_prime_prime_deriv(double w, void* data) {
  struct phi_prime_prime_params *params = (struct phi_prime_prime_params *)data;
  return pow(w, params->N_dot + params->k - 1) * exp(-w * params->H_dot) * 
    params->deriv_deriv_density(w, params->theta, params->deriv_idx_1, params->deriv_idx_2);
}

// The following phi functions are the only places that branch based on the
// frailty distribution. If a new distribution is implemented, it needs to go
// here, with the respective density and deriv density functions (or LT)

double phi(int k, int N_dot, double H_dot, double *theta, String frailty) {
  // Laplace transform integrals
  if (frailty == "gamma") {
    return lt_dgamma(N_dot + k - 1, H_dot, theta) * pow(-1, N_dot + k - 1);
  } else if (frailty == "posstab") {
    return lt_dposstab(N_dot + k - 1, H_dot, theta) * pow(-1, N_dot + k - 1);
  } else if (frailty == "pvf") {
    return lt_dpvf(N_dot + k - 1, H_dot, theta) * pow(-1, N_dot + k - 1);
  }
  // Numerical integration
  phi_params phi_p;
  if (frailty == "lognormal") {
    phi_p = (struct phi_params){k, N_dot, H_dot, theta, dlognormal};
  } else if (frailty == "invgauss") {
    phi_p = (struct phi_params){k, N_dot, H_dot, theta, dinvgauss};
  } else {
    throw std::range_error("Unsupported frailty distribution");
  }
  
  return integrate(&phi_deriv, &phi_p);
}

// phi using the derivative of the density wrt. parameter p[derive_idx]
double phi_prime(int k, int N_dot, double H_dot, double *theta, String frailty, int deriv_idx) {
  // Laplace transform integrals
  if (frailty == "gamma") {
    return deriv_lt_dgamma(N_dot + k - 1, H_dot, theta, deriv_idx) * pow(-1, N_dot + k - 1);
  } else if (frailty == "posstab") {
    return deriv_lt_dposstab(N_dot + k - 1, H_dot, theta, deriv_idx) * pow(-1, N_dot + k - 1);
  } else if (frailty == "pvf") {
    return deriv_lt_dpvf(N_dot + k - 1, H_dot, theta, deriv_idx) * pow(-1, N_dot + k - 1);
  }
  
  // Numerical integration
  phi_prime_params phi_p;
//   if (frailty == "gamma") {
//     phi_p = (struct phi_prime_params){k, N_dot, H_dot, theta, deriv_dgamma, deriv_idx};
//   } else 
  if (frailty == "lognormal") {
    phi_p = (struct phi_prime_params){k, N_dot, H_dot, theta, deriv_dlognormal, deriv_idx};
  } else if (frailty == "invgauss") {
    phi_p = (struct phi_prime_params){k, N_dot, H_dot, theta, deriv_dinvgauss, deriv_idx};
  } else {
    throw std::range_error("Unsupported frailty distribution");
  }
  
  return integrate(&phi_prime_deriv, &phi_p);
}

// phi using the derivative of the density wrt. parameter p[derive_idx]
double phi_prime_prime(int k, int N_dot, double H_dot, double *theta, String frailty, int deriv_idx_1, int deriv_idx_2) {
  // Laplace transform integrals
  if (frailty == "gamma") {
    return deriv_deriv_lt_dgamma(N_dot + k - 1, H_dot, theta, deriv_idx_1, deriv_idx_2) * pow(-1, N_dot + k - 1);
  } else if (frailty == "posstab") {
    throw std::range_error("Unsupported frailty distribution");
    // return deriv_lt_dposstab(N_dot + k - 1, H_dot, theta, deriv_idx) * pow(-1, N_dot + k - 1);
  } else if (frailty == "pvf") {
    throw std::range_error("Unsupported frailty distribution");
    // return deriv_lt_dpvf(N_dot + k - 1, H_dot, theta, deriv_idx) * pow(-1, N_dot + k - 1);
  }
  
  // Numerical integration
  phi_prime_prime_params phi_p;
  
  if (frailty == "lognormal") {
    throw std::range_error("Unsupported frailty distribution");
    // phi_p = (struct phi_prime_prime_params){k, N_dot, H_dot, theta, deriv_dlognormal, deriv_idx_1, deriv_idx_2};
  } else if (frailty == "invgauss") {
    throw std::range_error("Unsupported frailty distribution");
    // phi_p = (struct phi_prime_prime_params){k, N_dot, H_dot, theta, deriv_dinvgauss, deriv_idx_1, deriv_idx_2};
  } else {
    throw std::range_error("Unsupported frailty distribution");
  }
  
  return integrate(&phi_prime_prime_deriv, &phi_p);
}

double psi(int N_dot, double H_dot, double* theta, String frailty) {
  // Shortcut for gamma frailty
//   if (frailty == "gamma") {
//     return (N_dot + 1/theta[0])/(H_dot + 1/theta[0]);
//   }
  // return phi(2, N_dot, H_dot, theta, frailty)/phi(1, N_dot, H_dot, theta, frailty);
  double tmp1 = phi(1, N_dot, H_dot, theta, frailty);
  double tmp2 = phi(2, N_dot, H_dot, theta, frailty);
  // Rcout << tmp1 << " " << tmp2 << " " << N_dot << " " << H_dot << " " << theta << std::endl;
  return tmp2/tmp1;
}

// [[Rcpp::export]]
double phi_c(int k, int N_dot, double H_dot, double theta, String frailty) {
  return phi(k, N_dot, H_dot, &theta, frailty);
}

// [[Rcpp::export]]
double phi_prime_c(int k, int N_dot, double H_dot, double theta, String frailty, int deriv_idx) {
  return phi_prime(k, N_dot, H_dot, &theta, frailty, deriv_idx);
}

// [[Rcpp::export]]
Rcpp::List baseline_hazard_estimator(Rcpp::List X_, 
                                     Rcpp::List R_,
                                     NumericVector d_,
                                     Rcpp::List Y_, 
                                     Rcpp::List N_dot,
                                     NumericVector beta, 
                                     NumericVector theta, 
                                     String frailty) {
  int n_timesteps = d_.size();
  int n_clusters = X_.size();
  
  // H_ is a list of matrices
  Rcpp::List H_ = clone(Y_);
  
  // H_dot is a list of numeric vectors
  Rcpp::List H_dot = clone(N_dot);
  
  for (int i = 0; i < n_clusters; ++i) {
    // H_[i](NumericMatrix(n_members, n_timesteps));
    Rcpp::NumericMatrix H_i = H_(i);
    Rcpp::NumericVector H_dot_i = H_dot(i);
    for (int k = 0; k < n_timesteps; ++k) {
      for (int j = 0; j < H_i.nrow(); ++j) {
        H_i(j, k) = 0;
      }
      H_dot_i(k) = 0;
    }
  }
  // Lambda_hat is the baseline cumulative hazard estimate
  NumericVector Lambda_hat(n_timesteps);
  NumericVector delta_Lambda_hat(n_timesteps);
  
  double denom;
  // k starts from 1, not a typo
  for (int k = 1; k < n_timesteps; ++k) {
    denom = 0;
    for (int i = 0; i < n_clusters; ++i) {
      Rcpp::NumericMatrix X_i = X_(i);
      Rcpp::NumericMatrix Y_i = Y_(i);
      Rcpp::NumericVector N_dot_i = N_dot(i);
      Rcpp::NumericVector H_dot_i = H_dot(i);
      
      double tmp = 0;
      for (int j = 0; j < X_i.nrow(); ++j) {
        tmp += Y_i(j, k) * exp(sum(beta * X_i(j, _)));
      }
      
      // TODO: optimize by ignoring 0 terms, don't do the denom integrations
      denom += tmp * psi(N_dot_i(k - 1), 
                         H_dot_i(k - 1), 
                         theta.begin(), 
                         frailty);
    }
    
    delta_Lambda_hat(k) = d_(k)/denom;
    
    Lambda_hat(k) = Lambda_hat(k - 1) + delta_Lambda_hat(k);
    
    for (int i = 0; i < n_clusters; ++i) {
      Rcpp::NumericMatrix X_i = X_(i);
      Rcpp::NumericVector R_i = R_(i);
      Rcpp::NumericMatrix H_i = H_(i);
      Rcpp::NumericVector H_dot_i = H_dot(i);
      
      for (int j = 0; j < X_i.nrow(); ++j) {
        // R_ is the rank of failure as an R index
        int k_min = min(NumericVector::create(R_i(j) - 1, k));
        H_i(j, k) = Lambda_hat(k_min) * exp(sum(beta * X_i(j, _)));
      }
      
      H_dot_i(k) = sum(H_i( _ , k));
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("H_") = H_,
                            Rcpp::Named("H_dot") = H_dot,
                            Rcpp::Named("Lambda_hat") = Lambda_hat,
                            Rcpp::Named("lambda_hat") = delta_Lambda_hat);
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
    Rcpp::NumericVector R_i = R_(i); // Failure rank, whare T[r] = failure time
    Rcpp::NumericVector I_i = I_(i); // Failure indicator
    Rcpp::NumericVector N_dot_i = N_dot(i);
    Rcpp::NumericMatrix H_i = H_(i);
    Rcpp::NumericVector H_dot_i = H_dot(i);
    int tau_k = N_dot_i.size() - 1;
    
    double tmp = 0;
    for (int j = 0; j < X_i.nrow(); j++) {
      term1 += I_i[j] * X_i(j,beta_idx - 1);
      tmp += H_i(j, R_i[j] - 1) * X_i(j,beta_idx - 1);
    }
    term2 += tmp * psi(N_dot_i(tau_k), 
                       H_dot_i(tau_k), 
                       theta.begin(), 
                       frailty);
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
    Rcpp::NumericVector R_i = R_(i);
    Rcpp::NumericVector I_i = I_(i); // Failure indicator
    Rcpp::NumericVector N_dot_i = N_dot(i);
    Rcpp::NumericMatrix H_i = H_(i);
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
double loglikelihood(List X_,
                      List R_, 
                      List I_, 
                      List N_dot,
                      List H_dot,
                      NumericVector Lambda, 
                      NumericVector beta, 
                      NumericVector theta, 
                      String frailty) {
  
  int n_clusters = X_.size();
  
  double term1 = 0;
  double term2 = 0;
  for (int i = 0; i < n_clusters; ++i) {
    Rcpp::NumericMatrix X_i = X_(i);
    Rcpp::NumericVector I_i = I_(i); // Failure indicator
    Rcpp::NumericVector R_i = R_(i); // Failure rank, where T[r] = failure time
    Rcpp::NumericVector N_dot_i = N_dot(i);
    Rcpp::NumericVector H_dot_i = H_dot(i);
    int tau_k = N_dot_i.size() - 1;
    
    for (int j = 0; j < X_i.nrow(); j++) {
      // R_i[j] is an R index
      if (I_i(j) > 0) { 
        term1 += log(Lambda(R_i(j)-1) - Lambda(R_i(j)-2)) + sum(beta * X_i(j, _));
      }
    }
    
    term2 += log(phi(1, 
                     N_dot_i(tau_k), 
                     H_dot_i(tau_k), 
                     theta.begin(), 
                     frailty));
  }
  
  return term1 + term2;
}

// partial derivative of H_, H_dot, and Delta_Lambda_ wrt. beta
// Verified matched numerical results
// [[Rcpp::export]]
List dH_dbeta(NumericVector d_,
              List K_,
              List X_,
              List R_,
              List R_dot_, 
              List N_dot_,
              List H_, 
              List H_dot_,
              NumericVector Lambda,
              NumericVector Delta_Lambda,
              NumericVector beta,
              NumericVector theta,
              int beta_idx,
              String frailty) {
  beta_idx -= 1; // Convert R idx to C++ idx to avoid confusion
  int n_timesteps = d_.size();
  int n_clusters = X_.size();
  
  List dH_dbeta_ = clone(H_);
  List dH_dot_dbeta_ = clone(H_dot_);
  NumericVector dLambda_dbeta(n_timesteps);
  NumericVector dDelta_Lambda_dbeta(n_timesteps);
  
  for (int i = 0; i < n_clusters; ++i) {
    Rcpp::NumericMatrix dH_dbeta_i = dH_dbeta_(i);
    Rcpp::NumericVector dH_dot_dbeta_i = dH_dot_dbeta_(i);
    for (int k = 0; k < n_timesteps; ++k) {
      for (int j = 0; j < dH_dbeta_i.nrow(); ++j) {
        dH_dbeta_i(j, k) = 0;
      }
      dH_dot_dbeta_i(k) = 0;
    }
  }
  
  double factor1, factor2;
  
  // k starts from 1, not a typo
  for (int k = 1; k < n_timesteps; ++k) {
    factor1 = 0;
    factor2 = 0;
    
    for (int i = 0; i < n_clusters; ++i) {
      Rcpp::NumericMatrix X_i = X_(i);
      Rcpp::NumericMatrix R_i = R_(i);
      Rcpp::NumericVector R_dot_i = R_dot_(i);
      Rcpp::NumericVector N_dot_i = N_dot_(i);
      Rcpp::NumericVector H_dot_i = H_dot_(i);
      Rcpp::NumericVector dH_dot_dbeta_i = dH_dot_dbeta_(i);
      
      double factor2_term1 = dH_dot_dbeta_i(k - 1) * R_dot_i(k) *
        ((pow(phi(2,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty), 2)/
          pow(phi(1,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty), 2)) - 
         (phi(3,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty)/
          phi(1,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty)));
      
      double factor2_term2_tmp = 0;
      for (int j = 0; j < X_i.nrow(); ++j) {
        factor2_term2_tmp += R_i(j, k) * X_i(j, beta_idx);
      }
      
      double factor2_term2 = factor2_term2_tmp * 
        (phi(2,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty)/
         phi(1,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty));
      
      factor1 += R_dot_i(k) *
        phi(2,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty)/
        phi(1,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty);
      
      factor2 += factor2_term1 + factor2_term2;
    }
    dDelta_Lambda_dbeta(k) = (-d_(k)/pow(factor1, 2)) * factor2;
    dLambda_dbeta(k) = dLambda_dbeta(k - 1) + dDelta_Lambda_dbeta(k);
    
    for (int i = 0; i < n_clusters; ++i) {
      Rcpp::NumericMatrix X_i = X_(i);
      Rcpp::NumericVector K_i = K_(i);
      Rcpp::NumericMatrix dH_dbeta_i = dH_dbeta_(i);
      Rcpp::NumericVector dH_dot_dbeta_i = dH_dot_dbeta_(i);
      
      for (int j = 0; j < X_i.nrow(); ++j) {
        // K_ is the rank of failure as an R index
        int k_min = min(NumericVector::create(K_i(j) - 1, k));
        dH_dbeta_i(j, k) = dLambda_dbeta(k_min) * exp(sum(beta * X_i(j, _))) + 
          Lambda(k_min) * exp(sum(beta * X_i(j, _))) * X_i(j, beta_idx);
      }
      
      dH_dot_dbeta_i(k) = sum(dH_dbeta_i( _ , k));
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("dH_dbeta_") = dH_dbeta_,
                            Rcpp::Named("dH_dot_dbeta_") = dH_dot_dbeta_,
                            Rcpp::Named("dDelta_Lambda_dbeta") = dDelta_Lambda_dbeta,
                            Rcpp::Named("dLambda_dbeta") = dLambda_dbeta);
}

// partial derivative of H_, H_dot, and Delta_Lambda_ wrt. theta
// Verified matched numerical results
// [[Rcpp::export]]
List dH_dtheta(NumericVector d_,
              List K_,
              List X_,
              List R_,
              List R_dot_, 
              List N_dot_,
              List H_, 
              List H_dot_,
              NumericVector Lambda,
              NumericVector Delta_Lambda,
              NumericVector beta,
              NumericVector theta,
              int theta_idx,
              String frailty) {
  theta_idx -= 1; // Convert R idx to C++ idx to avoid confusion
  int n_timesteps = d_.size();
  int n_clusters = X_.size();
  
  List dH_dtheta_ = clone(H_);
  List dH_dot_dtheta_ = clone(H_dot_);
  NumericVector dLambda_dtheta(n_timesteps);
  NumericVector dDelta_Lambda_dtheta(n_timesteps);
  
  for (int i = 0; i < n_clusters; ++i) {
    Rcpp::NumericMatrix dH_dtheta_i = dH_dtheta_(i);
    Rcpp::NumericVector dH_dot_dtheta_i = dH_dot_dtheta_(i);
    for (int k = 0; k < n_timesteps; ++k) {
      for (int j = 0; j < dH_dtheta_i.nrow(); ++j) {
        dH_dtheta_i(j, k) = 0;
      }
      dH_dot_dtheta_i(k) = 0;
    }
  }
  
  double factor1, factor2;
  
  // k starts from 1, not a typo
  for (int k = 1; k < n_timesteps; ++k) {
    factor1 = 0;
    factor2 = 0;
    
    for (int i = 0; i < n_clusters; ++i) {
      Rcpp::NumericMatrix X_i = X_(i);
      Rcpp::NumericMatrix R_i = R_(i);
      Rcpp::NumericVector R_dot_i = R_dot_(i);
      Rcpp::NumericVector N_dot_i = N_dot_(i);
      Rcpp::NumericVector H_dot_i = H_dot_(i);
      Rcpp::NumericVector dH_dot_dtheta_i = dH_dot_dtheta_(i);
      
      double factor2_term1 = 
        phi_prime(2,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty, theta_idx)/
        phi(1,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty);
      
      double factor2_term2 = 
        phi(2,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty) *
        phi_prime(1,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty, theta_idx)/
        pow(phi(1,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty), 2);
      
      double factor2_term3_term1 =  
        pow(phi(2,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty), 2)/
        pow(phi(1,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty), 2);
      
      double factor2_term3_term2 =  
        phi(3,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty)/
        phi(1,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty);
      
      factor1 += R_dot_i(k) *
        phi(2,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty)/
        phi(1,N_dot_i(k - 1),H_dot_i(k - 1), theta.begin(), frailty);
      
      factor2 += R_dot_i(k) * (factor2_term1 - factor2_term2 + 
        dH_dot_dtheta_i(k - 1) * (factor2_term3_term1 - factor2_term3_term2));
    }
    dDelta_Lambda_dtheta(k) = (-d_(k)/pow(factor1, 2)) * factor2;
    dLambda_dtheta(k) = dLambda_dtheta(k - 1) + dDelta_Lambda_dtheta(k);
    
    for (int i = 0; i < n_clusters; ++i) {
      Rcpp::NumericMatrix X_i = X_(i);
      Rcpp::NumericVector K_i = K_(i);
      Rcpp::NumericMatrix dH_dtheta_i = dH_dtheta_(i);
      Rcpp::NumericVector dH_dot_dtheta_i = dH_dot_dtheta_(i);
      
      for (int j = 0; j < X_i.nrow(); ++j) {
        // K_ is the rank of failure as an R index
        int k_min = min(NumericVector::create(K_i(j) - 1, k));
        dH_dtheta_i(j, k) = dLambda_dtheta(k_min) * exp(sum(beta * X_i(j, _)));
      }
      
      dH_dot_dtheta_i(k) = sum(dH_dtheta_i( _ , k));
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("dH_dtheta_") = dH_dtheta_,
                            Rcpp::Named("dH_dot_dtheta_") = dH_dot_dtheta_,
                            Rcpp::Named("dDelta_Lambda_dtheta") = dDelta_Lambda_dtheta,
                            Rcpp::Named("dLambda_dtheta") = dLambda_dtheta);
}

// Eq. (34) in Gorfine (2008)
// [[Rcpp::export]]
double jacobian_beta_beta(Rcpp::List K_, 
                          Rcpp::List X_,
                          Rcpp::List N_dot,
                          Rcpp::List H_, 
                          Rcpp::List H_dot,
                          Rcpp::List dH_dbeta_, 
                          Rcpp::List dH_dot_dbeta_,
                          NumericVector beta, 
                          NumericVector theta,
                          int beta_idx_1, 
                          int beta_idx_2, 
                          String frailty) {
  beta_idx_1 -= 1; // Convert R idx to C++ idx to avoid confusion
  beta_idx_2 -= 1; 
  int n_clusters = X_.size();
  
  double out = 0;
  for (int i = 0; i < n_clusters; ++i) {
    Rcpp::NumericMatrix X_i = X_(i);
    Rcpp::NumericVector K_i = K_(i); // Failure rank, whare T[r] = failure time
    Rcpp::NumericVector N_dot_i = N_dot(i);
    Rcpp::NumericMatrix H_i = H_(i);
    Rcpp::NumericVector H_dot_i = H_dot(i);
    Rcpp::NumericMatrix dH_dbeta_i = dH_dbeta_(i);
    Rcpp::NumericVector dH_dot_dbeta_i = dH_dot_dbeta_(i);
    int tau_k = N_dot_i.size() - 1;
    
    // tmps for the inner loop sums
    double tmp_term1 = 0;
    double tmp_term2 = 0;
    
    for (int j = 0; j < X_i.nrow(); j++) {
      tmp_term1 += X_i(j, beta_idx_1) * dH_dbeta_i(j, K_i[j] - 1);
      tmp_term2 += H_i(j, K_i[j] - 1) * X_i(j, beta_idx_1) * dH_dot_dbeta_i(tau_k);
    }
    
    double term1 = tmp_term1 * 
      phi(2,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty)/
      phi(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty);
    
    double term2 = tmp_term2 * 
      ( (phi(3,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty)/
         phi(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty))
       -(pow(phi(2,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty),2)/
         pow(phi(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty),2)));
    
    out += term1 - term2;
  }
  out = -out/n_clusters;
  return out;
}

// Eq. (35) in Gorfine (2008)
// [[Rcpp::export]]
double jacobian_beta_theta(Rcpp::List K_, 
                          Rcpp::List X_,
                          Rcpp::List N_dot,
                          Rcpp::List H_, 
                          Rcpp::List H_dot,
                          Rcpp::List dH_dtheta_, 
                          Rcpp::List dH_dot_dtheta_,
                          NumericVector beta, 
                          NumericVector theta,
                          int beta_idx, 
                          int theta_idx, 
                          String frailty) {
  beta_idx -= 1; // Convert R idx to C++ idx to avoid confusion
  theta_idx -= 1;
  int n_clusters = X_.size();
  
  double out = 0;
  for (int i = 0; i < n_clusters; ++i) {
    Rcpp::NumericMatrix X_i = X_(i);
    Rcpp::NumericVector K_i = K_(i); // Failure rank
    Rcpp::NumericVector N_dot_i = N_dot(i);
    Rcpp::NumericMatrix H_i = H_(i);
    Rcpp::NumericVector H_dot_i = H_dot(i);
    Rcpp::NumericMatrix dH_dtheta_i = dH_dtheta_(i);
    Rcpp::NumericVector dH_dot_dtheta_i = dH_dot_dtheta_(i);
    int tau_k = N_dot_i.size() - 1;
    
    // tmps for the inner loop sums
    double tmp_term1 = 0;
    double tmp_term2 = 0;
    
    for (int j = 0; j < X_i.nrow(); j++) {
      tmp_term1 += X_i(j, beta_idx) * dH_dtheta_i(j, K_i[j] - 1);
      tmp_term2 += H_i(j, K_i[j] - 1) * X_i(j, beta_idx);
    }
    
    double term1 = tmp_term1 * 
      phi(2,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty)/
      phi(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty);
    
    double term2_term1 = 
      phi_prime(2,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty, theta_idx)/
      phi(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty);
    
    double term2_term2 = 
      phi(2,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty) *
      phi_prime(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty, theta_idx)/
      pow(phi(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty), 2);
    
    double term2_term3 = 
      ( (pow(phi(2,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty), 2)/
         pow(phi(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty), 2))
       -(phi(3,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty)/
         phi(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty))) *
         dH_dot_dtheta_i(tau_k);
    
    double term2 = (term2_term1 - term2_term2 + term2_term3) * tmp_term2;
    
    out += term1 + term2;
  }
  out = -out/n_clusters;
  return out;
}

// Eq. (36) in Gorfine (2008)
// [[Rcpp::export]]
double jacobian_theta_beta(Rcpp::List K_, 
                          Rcpp::List X_,
                          Rcpp::List N_dot,
                          Rcpp::List H_, 
                          Rcpp::List H_dot,
                          Rcpp::List dH_dbeta_, 
                          Rcpp::List dH_dot_dbeta_,
                          NumericVector beta, 
                          NumericVector theta,
                          int theta_idx, 
                          int beta_idx, 
                          String frailty) {
  theta_idx -= 1; // Convert R idx to C++ idx to avoid confusion
  beta_idx -= 1;
  int n_clusters = X_.size();
  
  double out = 0;
  for (int i = 0; i < n_clusters; ++i) {
    Rcpp::NumericMatrix X_i = X_(i);
    Rcpp::NumericVector K_i = K_(i); // Failure rank, whare T[r] = failure time
    Rcpp::NumericVector N_dot_i = N_dot(i);
    Rcpp::NumericMatrix H_i = H_(i);
    Rcpp::NumericVector H_dot_i = H_dot(i);
    Rcpp::NumericMatrix dH_dbeta_i = dH_dbeta_(i);
    Rcpp::NumericVector dH_dot_dbeta_i = dH_dot_dbeta_(i);
    int tau_k = N_dot_i.size() - 1;
    
    double term1 = 
      phi_prime(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty, theta_idx) *
      phi(2,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty)/
      pow(phi(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty), 2);
    
    double term2 = 
      phi_prime(2,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty, theta_idx)/
      phi(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty);
    
    out += (term1 - term2) * dH_dot_dbeta_i(tau_k);
  }
  out = out/n_clusters;
  return out;
}

// Eq. (37) in Gorfine (2008)
// [[Rcpp::export]]
double jacobian_theta_theta(Rcpp::List K_, 
                          Rcpp::List X_,
                          Rcpp::List N_dot,
                          Rcpp::List H_, 
                          Rcpp::List H_dot,
                          Rcpp::List dH_dtheta_, 
                          Rcpp::List dH_dot_dtheta_,
                          NumericVector beta, 
                          NumericVector theta,
                          int theta_idx_1, 
                          int theta_idx_2, 
                          String frailty) {
  theta_idx_1 -= 1; // Convert R idx to C++ idx to avoid confusion
  theta_idx_2 -= 1; 
  int n_clusters = X_.size();
  
  double out = 0;
  for (int i = 0; i < n_clusters; ++i) {
    Rcpp::NumericMatrix X_i = X_(i);
    Rcpp::NumericVector K_i = K_(i); // Failure rank, whare T[r] = failure time
    Rcpp::NumericVector N_dot_i = N_dot(i);
    Rcpp::NumericMatrix H_i = H_(i);
    Rcpp::NumericVector H_dot_i = H_dot(i);
    Rcpp::NumericMatrix dH_dtheta_i = dH_dtheta_(i);
    Rcpp::NumericVector dH_dot_dtheta_i = dH_dot_dtheta_(i);
    int tau_k = N_dot_i.size() - 1;
    
    double term1 =
      phi_prime_prime(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty, 
                      theta_idx_1, theta_idx_2)/
      phi(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty);
    
    double term2 = 
      pow(phi_prime(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty, theta_idx_1)/
          phi(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty), 2);
    
    // TODO: not sure about the theta indices
    double term3_factor1_term1 = 
      phi_prime(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty, theta_idx_2) *
      phi(2,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty)/
      pow(phi(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty), 2);
    
    double term3_factor1_term2 = 
      phi_prime(2,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty, theta_idx_2)/
      phi(1,N_dot_i(tau_k),H_dot_i(tau_k), theta.begin(), frailty);
    
    out += term1 - term2 + (term3_factor1_term1 - term3_factor1_term2) * dH_dot_dtheta_i(tau_k);
  }
  out = out/n_clusters;
  return out;
}


