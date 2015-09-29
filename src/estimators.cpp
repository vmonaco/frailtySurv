#include <Rcpp.h>
#include <cmath>
#include "cubature.h"
#include "distributions.h"

using namespace Rcpp;

// The following phi functions are the only places that branch based on the
// frailty distribution. If a new distribution is implemented, it needs to go
// here, with the Laplace transform functions

double phi(int k, int N_dot, double H_dot, double *theta, String frailty, 
           double abstol, double reltol, int maxit) {
  if (frailty == "gamma") {
    return lt_dgamma(N_dot + k - 1, H_dot, theta) * pow(-1.0, N_dot + k - 1);
  } else if (frailty == "pvf") {
    return lt_dpvf(N_dot + k - 1, H_dot, theta) * pow(-1.0, N_dot + k - 1);
  } else if (frailty == "lognormal") {
    return lt_dlognormal(N_dot + k - 1, H_dot, theta, abstol, reltol, maxit) * pow(-1.0, N_dot + k - 1);
  } else if (frailty == "invgauss") {
    return lt_dinvgauss(N_dot + k - 1, H_dot, theta, abstol, reltol, maxit) * pow(-1.0, N_dot + k - 1);
  } else {
    throw std::range_error("Unsupported frailty distribution");
  }
}

// phi using the derivative of the density wrt. parameter p[deriv_idx]
double phi_prime(int k, int N_dot, double H_dot, double *theta, String frailty, 
                 int deriv_idx, double abstol, double reltol, int maxit) {
  if (frailty == "gamma") {
    return deriv_lt_dgamma(N_dot + k - 1, H_dot, theta, deriv_idx) * pow(-1.0, N_dot + k - 1);
  } else if (frailty == "pvf") {
    return deriv_lt_dpvf(N_dot + k - 1, H_dot, theta, deriv_idx) * pow(-1.0, N_dot + k - 1);
  } else if (frailty == "lognormal") {
    return deriv_lt_dlognormal(N_dot + k - 1, H_dot, theta, deriv_idx, abstol, reltol, maxit) * pow(-1.0, N_dot + k - 1);
  } else if (frailty == "invgauss") {
    return deriv_lt_dinvgauss(N_dot + k - 1, H_dot, theta, deriv_idx, abstol, reltol, maxit) * pow(-1.0, N_dot + k - 1);
  } else {
    throw std::range_error("Unsupported frailty distribution");
  }
}

// phi using the derivative of the density wrt. parameter p[deriv_idx_1]
double phi_prime_prime(int k, int N_dot, double H_dot, double *theta, String frailty, 
                       int deriv_idx_1, int deriv_idx_2, double abstol, double reltol, int maxit) {
  if (frailty == "gamma") {
    return deriv_deriv_lt_dgamma(N_dot + k - 1, H_dot, theta, deriv_idx_1, deriv_idx_2) * pow(-1.0, N_dot + k - 1);
  } else if (frailty == "pvf") {
    return deriv_deriv_lt_dpvf(N_dot + k - 1, H_dot, theta, deriv_idx_1, deriv_idx_2) * pow(-1.0, N_dot + k - 1);
  } else if (frailty == "lognormal") {
    return deriv_deriv_lt_dlognormal(N_dot + k - 1, H_dot, theta, deriv_idx_1, deriv_idx_2, abstol, reltol, maxit) * pow(-1.0, N_dot + k - 1);
  } else if (frailty == "invgauss") {
    return deriv_deriv_lt_dinvgauss(N_dot + k - 1, H_dot, theta, deriv_idx_1, deriv_idx_2, abstol, reltol, maxit) * pow(-1.0, N_dot + k - 1);
  } else {
    throw std::range_error("Unsupported frailty distribution");
  }
}

// Mainly for testing, compare to the corresponding *_r functions
// [[Rcpp::export]]
double phi_c(int k, int N_dot, double H_dot, double theta, String frailty, 
             double abstol, double reltol, int maxit) {
  return phi(k, N_dot, H_dot, &theta, frailty, abstol, reltol, maxit);
}

// [[Rcpp::export]]
double phi_prime_c(int k, int N_dot, double H_dot, double theta, 
                   String frailty, int deriv_idx, 
                   double abstol, double reltol, int maxit) {
  return phi_prime(k, N_dot, H_dot, &theta, frailty, deriv_idx, 
                   abstol, reltol, maxit);
}

// [[Rcpp::export]]
double phi_prime_prime_c(int k, int N_dot, double H_dot, double theta, 
                         String frailty, int deriv_idx_1, int deriv_idx_2,
                         double abstol, double reltol, int maxit) {
  return phi_prime_prime(k, N_dot, H_dot, &theta, frailty, 
                         deriv_idx_1, deriv_idx_2,
                         abstol, reltol, maxit);
}

double psi(int N_dot, double H_dot, double* theta, String frailty, 
           double abstol, double reltol, int maxit) {
  // Shortcut for gamma frailty
//   if (frailty == "gamma") {
//     return (N_dot + 1/theta[0])/(H_dot + 1/theta[0]);
//   }
  return phi(2, N_dot, H_dot, theta, frailty, abstol, reltol, maxit)/
  phi(1, N_dot, H_dot, theta, frailty, abstol, reltol, maxit);
}

// Baseline hazard estimator, returns some other useful vars
// [[Rcpp::export]]
List bh(NumericVector d_,
        List R_star, 
        List K_,
        List Y_, 
        List N_,
        List N_dot,
        NumericVector beta, 
        NumericVector theta, 
        String frailty,
        NumericVector weights,
        double abstol, double reltol, int maxit) {
  
  int n_timesteps = d_.size();
  int n_clusters = R_star.size();
  
  // Lambda_hat is the baseline cumulative hazard estimate
  NumericVector Lambda(n_timesteps);
  NumericVector lambda(n_timesteps);
  
  // Y_ is a list of matrices
  List H_ = clone(Y_);
  
  // N_dot is a list of numeric vectors
  List H_dot = clone(N_dot);
  List psi_ = clone(N_dot);
  List phi_1_ = clone(N_dot);
  List phi_2_ = clone(N_dot);
  
  // Compute k=0
  for (int i = 0; i < n_clusters; ++i) {
    NumericMatrix H_i = H_(i);
    NumericVector H_dot_i = H_dot(i);
    NumericVector phi_1_i = phi_1_(i);
    NumericVector phi_2_i = phi_2_(i);
    NumericVector psi_i = psi_(i);
    
    for (int j = 0; j < H_i.nrow(); ++j) {
      H_i(j, 0) = 0;
    }
    H_dot_i(0) = 0;
    
    phi_1_i(0) = phi(1, 0, 0, theta.begin(), frailty, abstol, reltol, maxit);
    phi_2_i(0) = phi(2, 0, 0, theta.begin(), frailty, abstol, reltol, maxit);
    psi_i(0) = phi_2_i(0)/phi_1_i(0);
  }
  
  // k starts from 1, not a typo; values at k=0 are set above
  for (int k = 1; k < n_timesteps; ++k) {
    double denom = 0;
    double numer = 0;
    
    // Optimize by ignoring 0 terms, need to do the inner sums otherwise
    if (d_(k) > 0) {
      for (int i = 0; i < n_clusters; ++i) {
        NumericVector R_star_i = R_star(i);
        NumericMatrix Y_i = Y_(i);
        NumericMatrix N_i = N_(i);
        NumericVector N_dot_i = N_dot(i);
        NumericVector H_dot_i = H_dot(i);
        
        NumericVector psi_i = psi_(i);
        
        double denom_inner_sum = 0;
        double numer_inner_sum = 0;
        for (int j = 0; j < R_star_i.size(); ++j) {
          denom_inner_sum += Y_i(j, k) * R_star_i(j);
          numer_inner_sum += N_i(j, k) - N_i(j, k - 1);
        }
        
        denom += weights(i) * psi_i(k - 1) * denom_inner_sum;
        numer += weights(i) * numer_inner_sum;
      }
      
      lambda(k) = numer/denom;
    } else {
      lambda(k) = 0;
    }
    
    Lambda(k) = Lambda(k - 1) + lambda(k);
    
    for (int i = 0; i < n_clusters; ++i) {
      NumericVector R_star_i = R_star(i);
      NumericVector K_i = K_(i);
      NumericMatrix H_i = H_(i);
      NumericVector H_dot_i = H_dot(i);
      NumericVector N_dot_i = N_dot(i);
      
      NumericVector phi_1_i = phi_1_(i);
      NumericVector phi_2_i = phi_2_(i);
      NumericVector psi_i = psi_(i);
      
      
      // Updates for the next iteration
      if (d_(k) > 0) {
        for (int j = 0; j < R_star_i.size(); ++j) {
          // K_ is the rank of failure as an R index
          int k_min = min(NumericVector::create(K_i(j) - 1, k));
          H_i(j, k) = Lambda(k_min) * R_star_i(j);
        }
        
        H_dot_i(k) = sum(H_i( _ , k));
        phi_1_i(k) = phi(1, N_dot_i(k), H_dot_i(k), theta.begin(), frailty, abstol, reltol, maxit);
        phi_2_i(k) = phi(2, N_dot_i(k), H_dot_i(k), theta.begin(), frailty, abstol, reltol, maxit);
        psi_i(k) = phi_2_i(k)/phi_1_i(k);
      } else {
        H_i(_, k) = H_i(_, k - 1);
        H_dot_i(k) = H_dot_i(k - 1);
        phi_1_i(k) = phi_1_i(k - 1);
        phi_2_i(k) = phi_2_i(k - 1);
        psi_i(k) = psi_i(k - 1);
      }
    }
  }
  
  return List::create(Named("H_") = H_,
                      Named("H_dot") = H_dot,
                      Named("Lambda") = Lambda,
                      Named("lambda") = lambda,
                      Named("psi_") = psi_,
                      Named("phi_1_") = phi_1_,
                      Named("phi_2_") = phi_2_
  );
}

// [[Rcpp::export]]
NumericVector loglikelihood(List X_,
                     List K_, 
                     List I_, 
                     List phi_1_,
                     NumericVector lambda, 
                     NumericVector beta) {
  int k_tau = lambda.size() - 1;
  int n_clusters = X_.size();
  
  // double l = 0;
  NumericVector l(n_clusters);
  for (int i = 0; i < n_clusters; ++i) {
    
    NumericMatrix X_i = X_(i);
    NumericVector I_i = I_(i); // Failure indicator
    NumericVector K_i = K_(i); // R index
    NumericVector phi_1_i = phi_1_(i);
    
    double term1 = 0;
    for (int j = 0; j < X_i.nrow(); j++) {
      // K_i[j] is an R index
      if (I_i(j) > 0) { 
        term1 += log(lambda(K_i(j)-1)) + sum(beta * X_i(j, _));
      }
    }
    double term2 =  log(phi_1_i(k_tau));
    l(i) = term1 + term2;
  }
  
  return l;
}


// [[Rcpp::export]]
NumericVector xi_beta(List X_,
                      List I_, 
                      List H_,
                      List psi_,
                      int r) {
  r -= 1;  // Convert R idx to C++ idx to avoid confusion
  NumericVector tmp = psi_(0);
  int k_tau = tmp.size() - 1;
  int n_clusters = X_.size();
  NumericVector xi_r(n_clusters);
  
  for (int i = 0; i < n_clusters; ++i) {
    NumericMatrix X_i = X_(i);
    NumericMatrix H_i = H_(i);
    
    NumericVector I_i = I_(i);
    NumericVector psi_i = psi_(i);
    
    double term1 = 0;
    double term2_factor1 = 0;
    for (int j = 0; j < X_i.nrow(); j++) {
      term1 += I_i(j) * X_i(j, r);
      term2_factor1 += H_i(j, k_tau) * X_i(j, r);
    }
    
    double term2 = term2_factor1 * psi_i(k_tau);
    
    xi_r(i) = term1 - term2;
  }
  
  return xi_r;
}

// [[Rcpp::export]]
NumericVector xi_theta(List phi_1_,
                       List phi_prime_1_,
                       int r) {
  r -= 1;  // Convert R idx to C++ idx to avoid confusion
  NumericVector tmp = phi_1_(0);
  int k_tau = tmp.size() - 1;
  int n_clusters = phi_1_.size();
  NumericVector xi_r(n_clusters);
  
  for (int i = 0; i < n_clusters; ++i) {
    NumericVector phi_1_i = phi_1_(i);
    NumericVector phi_prime_1_i = phi_prime_1_(i);
    
    xi_r(i) = phi_prime_1_i(k_tau)/phi_1_i(k_tau);
  }
  
  return xi_r;
}

// [[Rcpp::export]]
List phi_k(int s,
           List N_dot_,
           List H_dot_,
           NumericVector theta, 
           String frailty,
           double abstol, double reltol, int maxit) {
  NumericVector tmp = N_dot_(0);
  int n_timesteps = tmp.size();
  int n_clusters = N_dot_.size();
  
  // N_dot is a list of numeric vectors
  List phi_k_ = clone(N_dot_);
  
  // Zero everything out
  for (int i = 0; i < n_clusters; ++i) {
    NumericVector N_dot_i = N_dot_(i);
    NumericVector H_dot_i = H_dot_(i);
    NumericVector phi_k_i = phi_k_(i);
    for (int k = 0; k < n_timesteps; ++k) {
        phi_k_i(k) = phi(s, N_dot_i(k), H_dot_i(k), theta.begin(), frailty, abstol, reltol, maxit);
    }
  }
  return phi_k_;
}

// [[Rcpp::export]]
List phi_prime_k(int s,
                 int theta_idx,
                 List N_dot_,
                 List H_dot_,
                 NumericVector theta, 
                 String frailty,
                 double abstol, double reltol, int maxit) {
  theta_idx -= 1;
  NumericVector tmp = N_dot_(0);
  int n_timesteps = tmp.size();
  int n_clusters = N_dot_.size();
  
  // N_dot is a list of numeric vectors
  List phi_prime_k_ = clone(N_dot_);
  
  // Zero everything out
  for (int i = 0; i < n_clusters; ++i) {
    NumericVector N_dot_i = N_dot_(i);
    NumericVector H_dot_i = H_dot_(i);
    NumericVector phi_prime_k_i = phi_prime_k_(i);
    for (int k = 0; k < n_timesteps; ++k) {
      phi_prime_k_i(k) = phi_prime(s, N_dot_i(k), H_dot_i(k),
                    theta.begin(), frailty, theta_idx,
                    abstol, reltol, maxit);
    }
  }
  return phi_prime_k_;
}

// [[Rcpp::export]]
List phi_prime_prime_k(int s,
                       int theta_idx_1,
                       int theta_idx_2,
                       List N_dot_,
                       List H_dot_,
                       NumericVector theta, 
                       String frailty,
                       int kstart,
                       double abstol, double reltol, int maxit) {
  theta_idx_1 -= 1;
  theta_idx_2 -= 1;
  kstart -= 1;
  NumericVector tmp = N_dot_(0);
  int n_timesteps = tmp.size();
  int n_clusters = N_dot_.size();
  
  // N_dot is a list of numeric vectors
  List phi_prime_prime_k_ = clone(N_dot_);
  
  for (int i = 0; i < n_clusters; ++i) {
    NumericVector N_dot_i = N_dot_(i);
    NumericVector H_dot_i = H_dot_(i);
    NumericVector phi_prime_prime_k_i = phi_prime_prime_k_(i);
    for (int k = kstart; k < n_timesteps; ++k) {
      phi_prime_prime_k_i(k) = phi_prime_prime(s, N_dot_i(k), H_dot_i(k), 
                          theta.begin(), frailty, theta_idx_1, theta_idx_2,
                          abstol, reltol, maxit);
    }
  }
  return phi_prime_prime_k_;
}

// partial derivative of H_, H_dot_, and Lambda wrt. beta
// Verified these match numerical results
// [[Rcpp::export]]
List dH_dbeta(int s,
              NumericVector d_,
              List X_,
              List K_,
              List R_,
              List R_dot_,
              List R_star,
              List phi_1_,
              List phi_2_,
              List phi_3_,
              NumericVector Lambda,
              NumericVector lambda,
              NumericVector beta,
              NumericVector theta,
              String frailty) {
  s -= 1; // Convert R idx to C++ idx to avoid confusion
  int n_timesteps = d_.size();
  int n_clusters = X_.size();
  
  List dH_dbeta_ = clone(R_);
  List dH_dot_dbeta_ = clone(R_dot_);
  NumericVector dLambda_dbeta(n_timesteps);
  NumericVector dlambda_dbeta(n_timesteps);
  
  for (int i = 0; i < n_clusters; ++i) {
    NumericMatrix dH_dbeta_i = dH_dbeta_(i);
    NumericVector dH_dot_dbeta_i = dH_dot_dbeta_(i);
    for (int k = 0; k < n_timesteps; ++k) {
      for (int j = 0; j < dH_dbeta_i.nrow(); ++j) {
        dH_dbeta_i(j, k) = 0;
      }
      dH_dot_dbeta_i(k) = 0;
    }
  }
  
  // k starts from 1, not a typo
  for (int k = 1; k < n_timesteps; ++k) {
    double factor1 = 0;
    double factor2 = 0;
    
    for (int i = 0; i < n_clusters; ++i) {
      NumericMatrix X_i = X_(i);
      NumericMatrix R_i = R_(i);
      NumericVector R_dot_i = R_dot_(i);
      NumericVector phi_1_i = phi_1_(i);
      NumericVector phi_2_i = phi_2_(i);
      NumericVector phi_3_i = phi_3_(i);
      NumericVector dH_dot_dbeta_i = dH_dot_dbeta_(i);
      
      double factor2_term1 = dH_dot_dbeta_i(k - 1) * R_dot_i(k) *
        ((pow(phi_2_i(k - 1), 2)/
          pow(phi_1_i(k - 1), 2)) - 
         (phi_3_i(k - 1)/
          phi_1_i(k - 1)));
      
      double factor2_term2_tmp = 0;
      for (int j = 0; j < X_i.nrow(); ++j) {
        factor2_term2_tmp += R_i(j, k) * X_i(j, s);
      }
      
      double factor2_term2 = factor2_term2_tmp * (phi_2_i(k - 1)/phi_1_i(k - 1));
      
      factor1 += R_dot_i(k) * (phi_2_i(k - 1)/phi_1_i(k - 1));
      factor2 += factor2_term1 + factor2_term2;
    }
    dlambda_dbeta(k) = (-d_(k)/pow(factor1, 2)) * factor2;
    dLambda_dbeta(k) = dLambda_dbeta(k - 1) + dlambda_dbeta(k);
    
    for (int i = 0; i < n_clusters; ++i) {
      NumericMatrix X_i = X_(i);
      NumericVector R_star_i = R_star(i);
      NumericVector K_i = K_(i);
      NumericMatrix dH_dbeta_i = dH_dbeta_(i);
      NumericVector dH_dot_dbeta_i = dH_dot_dbeta_(i);
      
      for (int j = 0; j < X_i.nrow(); ++j) {
        // K_ is the rank of failure as an R index
        int k_min = min(NumericVector::create(K_i(j) - 1, k));
        dH_dbeta_i(j, k) = dLambda_dbeta(k_min) * R_star_i(j) + 
          Lambda(k_min) * R_star_i(j) * X_i(j, s);
      }
      
      dH_dot_dbeta_i(k) = sum(dH_dbeta_i( _ , k));
    }
  }
  
  return List::create(Named("dH_dbeta_") = dH_dbeta_,
                            Named("dH_dot_dbeta_") = dH_dot_dbeta_,
                            Named("dlambda_dbeta") = dlambda_dbeta,
                            Named("dLambda_dbeta") = dLambda_dbeta);
}

// partial derivative of H_, H_dot, and lambda_ wrt. theta
// Verified matched numerical results
// [[Rcpp::export]]
List dH_dtheta(NumericVector d_,
               List X_,
               List K_,
               List R_,
               List R_dot_,
               List R_star,
               List phi_1_,
               List phi_2_,
               List phi_3_,
               List phi_prime_1_,
               List phi_prime_2_,
               NumericVector Lambda,
               NumericVector lambda,
               NumericVector beta) {
  int n_timesteps = d_.size();
  int n_clusters = X_.size();
  
  List dH_dtheta_ = clone(R_);
  List dH_dot_dtheta_ = clone(R_dot_);
  NumericVector dLambda_dtheta(n_timesteps);
  NumericVector dlambda_dtheta(n_timesteps);
  
  for (int i = 0; i < n_clusters; ++i) {
    NumericMatrix dH_dtheta_i = dH_dtheta_(i);
    NumericVector dH_dot_dtheta_i = dH_dot_dtheta_(i);
    for (int k = 0; k < n_timesteps; ++k) {
      for (int j = 0; j < dH_dtheta_i.nrow(); ++j) {
        dH_dtheta_i(j, k) = 0;
      }
      dH_dot_dtheta_i(k) = 0;
    }
  }
  
  // k starts from 1, not a typo
  for (int k = 1; k < n_timesteps; ++k) {
    double factor1 = 0;
    double factor2 = 0;
    
    if (d_(k) > 0) {
      for (int i = 0; i < n_clusters; ++i) {
        NumericVector R_dot_i = R_dot_(i);
        
        NumericVector phi_1_i = phi_1_(i);
        NumericVector phi_2_i = phi_2_(i);
        NumericVector phi_3_i = phi_3_(i);
        NumericVector phi_prime_1_i = phi_prime_1_(i);
        NumericVector phi_prime_2_i = phi_prime_2_(i);
        
        NumericVector dH_dot_dtheta_i = dH_dot_dtheta_(i);
        
        double factor2_term1 = phi_prime_2_i(k - 1)/phi_1_i(k - 1);
        double factor2_term2 = phi_2_i(k - 1) * phi_prime_1_i(k - 1)/
        pow(phi_1_i(k - 1), 2);
        double factor2_term3_term1 = pow(phi_2_i(k - 1), 2)/pow(phi_1_i(k - 1), 2);
        double factor2_term3_term2 = phi_3_i(k - 1)/phi_1_i(k - 1);
        
        factor1 += R_dot_i(k) * phi_2_i(k - 1)/phi_1_i(k - 1);
        factor2 += R_dot_i(k) * (factor2_term1 - factor2_term2 + 
          dH_dot_dtheta_i(k - 1) * (factor2_term3_term1 - factor2_term3_term2));
      }
      dlambda_dtheta(k) = (-d_(k)/pow(factor1, 2)) * factor2;
    } else {
      dlambda_dtheta(k) = 0;
    }
    
    dLambda_dtheta(k) = dLambda_dtheta(k - 1) + dlambda_dtheta(k);
    
    for (int i = 0; i < n_clusters; ++i) {
      NumericMatrix X_i = X_(i);
      NumericVector R_star_i = R_star(i);
      NumericVector K_i = K_(i);
      NumericMatrix dH_dtheta_i = dH_dtheta_(i);
      NumericVector dH_dot_dtheta_i = dH_dot_dtheta_(i);
      
      for (int j = 0; j < X_i.nrow(); ++j) {
        // K_ is the rank of failure as an R index
        int k_min = min(NumericVector::create(K_i(j) - 1, k));
        dH_dtheta_i(j, k) = dLambda_dtheta(k_min) * R_star_i(j);
      }
      dH_dot_dtheta_i(k) = sum(dH_dtheta_i( _ , k));
    }
  }
  
  return List::create(Named("dH_dtheta_") = dH_dtheta_,
                            Named("dH_dot_dtheta_") = dH_dot_dtheta_,
                            Named("dlambda_dtheta") = dlambda_dtheta,
                            Named("dLambda_dtheta") = dLambda_dtheta);
}

// Eq. (34) in Gorfine (2008)
// [[Rcpp::export]]
double jacobian_beta_beta(int l,
                          List X_,
                          List K_, 
                          List H_, 
                          List phi_1_,
                          List phi_2_,
                          List phi_3_,
                          List dH_dbeta_, 
                          List dH_dot_dbeta_) {
  l -= 1; // Convert R idx to C++ idx to avoid confusion
  NumericVector tmp = dH_dot_dbeta_(0);
  int k_tau = tmp.size() - 1;
  int n_clusters = X_.size();
  
  double out = 0;
  for (int i = 0; i < n_clusters; ++i) {
    NumericVector K_i = K_(i); // Failure rank, whare T[k] = failure time
    NumericMatrix X_i = X_(i);
    NumericMatrix H_i = H_(i);
    
    NumericVector phi_1_i = phi_1_(i);
    NumericVector phi_2_i = phi_2_(i);
    NumericVector phi_3_i = phi_3_(i);
    
    NumericMatrix dH_dbeta_i = dH_dbeta_(i);
    NumericVector dH_dot_dbeta_i = dH_dot_dbeta_(i);
    
    // tmps for the inner loop sums
    double tmp_term1 = 0;
    double tmp_term2 = 0;
    
    for (int j = 0; j < X_i.nrow(); j++) {
      tmp_term1 += X_i(j, l) * dH_dbeta_i(j, K_i(j) - 1);
      tmp_term2 += H_i(j, K_i(j) - 1) * X_i(j, l) * dH_dot_dbeta_i(k_tau);
    }
    
    double term1 = tmp_term1 * phi_2_i(k_tau)/phi_1_i(k_tau);
    
    double term2 = tmp_term2 * 
      ( (phi_3_i(k_tau)/
         phi_1_i(k_tau))
       -(pow(phi_2_i(k_tau),2)/
         pow(phi_1_i(k_tau),2)));
    
    out += term1 - term2;
  }
  out = -out/n_clusters;
  return out;
}

// Eq. (35) in Gorfine (2008)
// [[Rcpp::export]]
double jacobian_beta_theta(int l,
                           List X_,
                           List K_, 
                           List H_, 
                           List phi_1_,
                           List phi_2_,
                           List phi_3_,
                           List phi_prime_1_,
                           List phi_prime_2_,
                           List dH_dtheta_, 
                           List dH_dot_dtheta_) {
  l -= 1; // Convert R idx to C++ idx to avoid confusion
  NumericVector tmp = dH_dot_dtheta_(0);
  int k_tau = tmp.size() - 1;
  int n_clusters = X_.size();
  
  double out = 0;
  for (int i = 0; i < n_clusters; ++i) {
    NumericVector K_i = K_(i); // Failure rank
    NumericMatrix X_i = X_(i);
    NumericMatrix H_i = H_(i);
    
    NumericVector phi_1_i = phi_1_(i);
    NumericVector phi_2_i = phi_2_(i);
    NumericVector phi_3_i = phi_3_(i);
    NumericVector phi_prime_1_i = phi_prime_1_(i);
    NumericVector phi_prime_2_i = phi_prime_2_(i);
    
    NumericMatrix dH_dtheta_i = dH_dtheta_(i);
    NumericVector dH_dot_dtheta_i = dH_dot_dtheta_(i);
    
    // tmps for the inner loop sums
    double tmp_term1 = 0;
    double tmp_term2 = 0;
    
    for (int j = 0; j < X_i.nrow(); j++) {
      tmp_term1 += X_i(j, l) * dH_dtheta_i(j, K_i(j) - 1);
      tmp_term2 += H_i(j, K_i(j) - 1) * X_i(j, l);
    }
    
    double term1 = tmp_term1 * phi_2_i(k_tau)/phi_1_i(k_tau);
    
    double term2_term1 = phi_prime_2_i(k_tau)/phi_1_i(k_tau);
    
    double term2_term2 = phi_2_i(k_tau) * phi_prime_1_i(k_tau)/
                                          pow(phi_1_i(k_tau), 2);
    
    double term2_term3 = dH_dot_dtheta_i(k_tau) *
      ( (pow(phi_2_i(k_tau), 2)/
         pow(phi_1_i(k_tau), 2))
       -(phi_3_i(k_tau)/
         phi_1_i(k_tau)));
    
    double term2 = (term2_term1 - term2_term2 + term2_term3) * tmp_term2;
    out += term1 + term2;
  }
  out = -out/n_clusters;
  return out;
}

// Eq. (36) in Gorfine (2008)
// [[Rcpp::export]]
double jacobian_theta_beta(List phi_1_,
                           List phi_2_,
                           List phi_prime_1_,
                           List phi_prime_2_,
                           List dH_dot_dbeta_) {
  NumericVector tmp = dH_dot_dbeta_(0);
  int k_tau = tmp.size() - 1;
  int n_clusters = dH_dot_dbeta_.size();
  
  double out = 0;
  for (int i = 0; i < n_clusters; ++i) {
    NumericVector phi_1_i = phi_1_(i);
    NumericVector phi_2_i = phi_2_(i);
    NumericVector phi_prime_1_i = phi_prime_1_(i);
    NumericVector phi_prime_2_i = phi_prime_2_(i);
    
    NumericVector dH_dot_dbeta_i = dH_dot_dbeta_(i);
    
    double term1 = phi_prime_1_i(k_tau) * phi_2_i(k_tau)/pow(phi_1_i(k_tau), 2);
    double term2 = phi_prime_2_i(k_tau)/phi_1_i(k_tau);
    
    out += (term1 - term2) * dH_dot_dbeta_i(k_tau);
  }
  out = out/n_clusters;
  return out;
}

// Eq. (37) in Gorfine (2008)
// [[Rcpp::export]]
double jacobian_theta_theta(List phi_1_,
                            List phi_2_,
                            List phi_prime_1_,
                            List phi_prime_2_,
                            List phi_prime_prime_1_,
                            List dH_dot_dtheta_) {
  NumericVector tmp = dH_dot_dtheta_(0);
  int k_tau = tmp.size() - 1;
  int n_clusters = dH_dot_dtheta_.size();
  
  double out = 0;
  for (int i = 0; i < n_clusters; ++i) {
    NumericVector phi_1_i = phi_1_(i);
    NumericVector phi_2_i = phi_2_(i);
    NumericVector phi_prime_1_i = phi_prime_1_(i);
    NumericVector phi_prime_2_i = phi_prime_2_(i);
    NumericVector phi_prime_prime_1_i = phi_prime_prime_1_(i);
    
    NumericVector dH_dot_dtheta_i = dH_dot_dtheta_(i);
    
    double term1 =phi_prime_prime_1_i(k_tau)/phi_1_i(k_tau);
    double term2 = pow(phi_prime_1_i(k_tau)/phi_1_i(k_tau), 2);
    double term3_factor1_term1 = phi_prime_1_i(k_tau) * phi_2_i(k_tau)/pow(phi_1_i(k_tau), 2);
    double term3_factor1_term2 = phi_prime_2_i(k_tau)/phi_1_i(k_tau);
    
    out += term1 - term2 + (term3_factor1_term1 - term3_factor1_term2) * dH_dot_dtheta_i(k_tau);
  }
  out = out/n_clusters;
  return out;
}

// [[Rcpp::export]]
List Q_beta(List X_,
            List K_,
            List H_, 
            List R_star,
            List phi_1_,
            List phi_2_,
            List phi_3_,
            int r) {
  r -= 1; // Convert R idx to C++ idx to avoid confusion
  NumericVector tmp = phi_1_(0);
  int n_timesteps = tmp.size();
  int n_clusters = X_.size();
  int k_tau = n_timesteps - 1;
  
  List Q_beta_ = clone(H_);
  
  for (int i = 0; i < n_clusters; ++i) {
    // Output
    NumericMatrix Q_beta_i = Q_beta_(i);
    
    // Matrices
    NumericMatrix X_i = X_(i);
    NumericMatrix H_i = H_(i);
    NumericVector K_i = K_(i);
    
    // Vectors
    NumericVector R_star_i = R_star(i);
    NumericVector phi_1_i = phi_1_(i);
    NumericVector phi_2_i = phi_2_(i);
    NumericVector phi_3_i = phi_3_(i);
  
    // TODO: j is reassigned in eq. on page 24?
    for (int j = 0; j < Q_beta_i.nrow(); ++j) {
      for (int k = 0; k < n_timesteps; ++k) {
        double inner_sum = 0;
        for (int l = 0; l < Q_beta_i.nrow(); ++l) {
          inner_sum += H_i(l, k) * X_i(l, r);
          // inner_sum += H_i(l, K_i(l)-1) * X_i(l, r);
        }
        
        double term1 = X_i(j, r) * (phi_2_i(k_tau)/phi_1_i(k_tau));
        double term2 = inner_sum * (phi_3_i(k_tau)/phi_1_i(k_tau));
        double term3 = inner_sum * (pow(phi_2_i(k_tau), 2)/
                                    pow(phi_1_i(k_tau), 2));
          
        Q_beta_i(j, k) = -R_star_i(j) * (term1 - term2 + term3);
      }
    }
  }
  
  return Q_beta_;
}

// [[Rcpp::export]]
List Q_theta(List H_,
             List R_star,
             List phi_1_,
             List phi_2_,
             List phi_prime_1_,
             List phi_prime_2_) {
  NumericVector tmp = phi_1_(0);
  int n_timesteps = tmp.size();
  int n_clusters = phi_1_.size();
  int k_tau = n_timesteps - 1;
  
  List Q_theta_ = clone(H_);
  
  for (int i = 0; i < n_clusters; ++i) {
    // Output
    NumericMatrix Q_theta_i = Q_theta_(i);
    
    // Vectors
    NumericVector R_star_i = R_star(i);
    NumericVector phi_1_i = phi_1_(i);
    NumericVector phi_2_i = phi_2_(i);
    NumericVector phi_prime_1_i = phi_prime_1_(i);
    NumericVector phi_prime_2_i = phi_prime_2_(i);
    
    double term1 = phi_2_i(k_tau) * phi_prime_1_i(k_tau)/pow(phi_1_i(k_tau), 2);
    double term2 = phi_prime_2_i(k_tau)/phi_1_i(k_tau);
    
    for (int j = 0; j < Q_theta_i.nrow(); ++j) {
      for (int k = 0; k < n_timesteps; ++k) {
        Q_theta_i(j, k) = R_star_i(j) * (term1 - term2);
      }
    }
  }
  
  return Q_theta_;
}

// [[Rcpp::export]]
NumericVector Ycal(List X_,
                   List R_star,
                   List Y_,
                   List psi_,
                   NumericVector beta) {
  NumericVector tmp = psi_(0);
  int n_timesteps = tmp.size();
  int n_clusters = X_.size();
  int k_tau = n_timesteps - 1;
  
  NumericVector Ycal_(n_timesteps);
  NumericVector Upsilon_(n_timesteps);
  
  for (int s = 0; s < n_timesteps; ++s) {
    double outer_sum = 0;
    for (int i = 0; i < n_clusters; ++i) {
      // Matrices
      NumericMatrix X_i = X_(i);
      NumericMatrix Y_i = Y_(i);
      
      // Vectors
      NumericVector R_star_i = R_star(i);
      NumericVector psi_i = psi_(i);
      
      double inner_sum = 0;
      for (int j = 0; j < X_i.nrow(); ++j) {
        inner_sum += Y_i(j, s) * R_star_i(j);
      }
      
      outer_sum += psi_i(s) * inner_sum;
    }
    Ycal_(s) = outer_sum/n_clusters;
  }
  
  return Ycal_;
}

// [[Rcpp::export]]
List eta(List phi_1_,
         List phi_2_,
         List phi_3_) {
  NumericVector tmp = phi_1_(0);
  int n_timesteps = tmp.size();
  int n_clusters = phi_1_.size();
  int k_tau = n_timesteps - 1;
  
  List eta_ = clone(phi_1_);
  
  for (int i = 0; i < n_clusters; ++i) {
    // Output
    NumericVector eta_i = eta_(i);
    
    // Vectors
    NumericVector phi_1_i = phi_1_(i);
    NumericVector phi_2_i = phi_2_(i);
    NumericVector phi_3_i = phi_3_(i);
      
    for (int s = 0; s < n_timesteps; ++s) {
      eta_i(s) = (phi_3_i(s)/phi_1_i(s)) - pow(phi_2_i(s)/phi_1_i(s), 2);
    }
  }
  
  return eta_;
}

// [[Rcpp::export]]
NumericVector Upsilon(List X_,
                      List R_star,
                      List K_,
                      List R_dot_,
                      List eta_,
                      NumericVector Ycal_,
                      NumericVector beta) {
  NumericVector tmp = R_dot_(0);
  int n_timesteps = tmp.size();
  int n_clusters = X_.size();
  
  NumericVector Upsilon_(n_timesteps);
  
  for (int s = 0; s < n_timesteps; ++s) {
    
    double factor3 = 0;
    for (int i = 0; i < n_clusters; ++i) {
      // Matrices
      NumericMatrix X_i = X_(i);
      
      // Vectors
      NumericVector R_star_i = R_star(i);
      NumericVector K_i = K_(i);
      NumericVector R_dot_i = R_dot_(i);
      NumericVector eta_i = eta_(i);
      
      for (int j = 0; j < X_i.nrow(); ++j) {
        if ((K_i(j)-1) > s) {
          factor3 += R_dot_i(s) * eta_i(s) * R_star_i(j);
        }
      }
    }
    
    Upsilon_(s) = pow((double)n_clusters, -2) * pow(Ycal_(s), -2) * factor3;
  }
  
  return Upsilon_;
}

// [[Rcpp::export]]
List Omega(List X_,
           List R_star,
           List N_,
           List R_dot_,
           List eta_,
           NumericVector Ycal_,
           NumericVector beta) {
  NumericVector tmp = R_dot_(0);
  int n_timesteps = tmp.size();
  int n_clusters = N_.size();
  List Omega_ = clone(N_);
  
  // Precompute dN_ij sum over i,j
  NumericVector dN_sum(n_timesteps);
  for (int u = 1; u < n_timesteps; ++u) {
    for (int i = 0; i < n_clusters; ++i) {
      NumericMatrix N_i = N_(i);
      for (int j = 0; j < N_i.nrow(); ++j) {
        dN_sum(u) += N_i(j, u) - N_i(j, u - 1);
      }
    }
  }
  
  double integration;
  for (int i = 0; i < n_clusters; ++i) {
    NumericMatrix X_i = X_(i);
    NumericMatrix N_i = N_(i);
    NumericMatrix Omega_i = Omega_(i);
    
    NumericVector R_star_i = R_star(i);
    NumericVector R_dot_i = R_dot_(i);
    NumericVector eta_i = eta_(i);
    
    for (int j = 0; j < Omega_i.nrow(); ++j) {
      
      for (int u = 1; u < n_timesteps; ++u) {
        Omega_i(j, u) = pow(Ycal_(u), -2) * R_dot_i(u) * eta_i(u) * 
          R_star_i(j) * dN_sum(u);
      }
    }
  }
  
  return Omega_;
}

// [[Rcpp::export]]
NumericVector p_hat(List I_,
                    NumericVector Upsilon_,
                    List Omega_,
                    List N_tilde_) {
  int n_timesteps = Upsilon_.size();
  int n_clusters = N_tilde_.size();
  
  NumericVector p_hat_(n_timesteps);
  
  for (int t = 0; t < n_timesteps; ++t) {
    double product = 1;
    
    for (int s = 1; s <= t; ++s) {
      double inner_sum = 0;
      
      for (int i = 0; i < n_clusters; ++i) {
        NumericVector I_i = I_(i);
        NumericMatrix N_tilde_i = N_tilde_(i);
        NumericMatrix Omega_i = Omega_(i);
        
        for (int j = 0; j < N_tilde_i.nrow(); ++j) {
          
          double omega_integral = 0;
          for (int u = s; u <= t; ++u) omega_integral += Omega_i(j, u);
          omega_integral = pow((double)n_clusters, -2) * omega_integral;
          
          inner_sum += (I_i(j) * Upsilon_(s) + omega_integral) * 
            (N_tilde_i(j, s) - N_tilde_i(j, s - 1));
        }
      }
      
      product *= 1 + inner_sum;
    }
    
    p_hat_(t) = product;
  }
  
  
  return p_hat_;
}

// [[Rcpp::export]]
NumericVector pi_r(List Q_,
                   List N_tilde_,
                   NumericVector p_hat) {
  
  int n_timesteps = p_hat.size();
  int k_tau = n_timesteps - 1;
  int n_clusters = N_tilde_.size();
  
  NumericVector pi_r_(n_timesteps);
  
  double integration, inner_sum;
  for (int s = 0; s < n_timesteps; ++s) {
    integration = 0;
    
    for (int t = s+1; t < n_timesteps; ++t) {
      inner_sum = 0;
      
      for (int i = 0; i < n_clusters; ++i) {
        NumericMatrix Q_i = Q_(i);
        NumericMatrix N_tilde_i = N_tilde_(i);
        
        for (int j = 0; j < N_tilde_i.nrow(); ++j) {
          inner_sum += Q_i(j, t) * (N_tilde_i(j, t) - N_tilde_i(j, t - 1));
        }
      }
      integration += inner_sum/p_hat(t);
    }
    
    pi_r_(s) = integration/n_clusters;
  }
  
  return pi_r_;
}

// [[Rcpp::export]]
double G_rl(NumericVector pi_r,
            NumericVector pi_l,
            NumericVector p_hat,
            NumericVector Ycal_,
            List N_) {
  
  int n_timesteps = p_hat.size();
  int k_tau = n_timesteps - 1;
  int n_clusters = N_.size();
  
  double integration = 0;
  for (int s = 1; s < n_timesteps; ++s) {
    
    double inner_sum = 0;
    for (int i = 0; i < n_clusters; ++i) {
      NumericMatrix N_i = N_(i);
      for (int j = 0; j < N_i.nrow(); ++j) {
        inner_sum += (N_i(j, s) - N_i(j, s - 1));
      }
    }
    
    integration += pi_r(s) * pi_l(s) * pow(p_hat(s - 1), 2) * inner_sum/pow(Ycal_(s), 2);
  }
  
  return integration/n_clusters;
}

// [[Rcpp::export]]
List M_hat(List X_,
           List R_star,
           List N_,
           List Y_,
           List psi_,
           NumericVector beta,
           NumericVector Lambda) {
  
  int n_timesteps = Lambda.size();
  int n_clusters = X_.size();
  
  List M_hat_ = clone(N_);
  
  for (int i = 0; i < n_clusters; ++i) {
    NumericMatrix X_i = X_(i);
    NumericMatrix N_i = N_(i);
    NumericMatrix Y_i = Y_(i);
    NumericMatrix M_hat_i = M_hat_(i);
    NumericVector R_star_i = R_star(i);
    NumericVector psi_i = psi_(i);
    
    for (int j = 0; j < X_i.nrow(); ++j) {
      for (int t = 0; t < n_timesteps; ++t) {
        double integration = 0;
        
        for (int u = 1; u <= t; ++u) {
          integration += R_star_i(j) * Y_i(j, u) * psi_i(u - 1) * 
            (Lambda(u) - Lambda(u - 1));
        }
        
        M_hat_i(j, t) = N_i(j, t) - integration;
      }
    }
  }
  
  return M_hat_;
}

// [[Rcpp::export]]
NumericMatrix u_star(List pi_,
                     NumericVector p_hat,
                     NumericVector Ycal_,
                     List M_hat_) {
  int n_timesteps = p_hat.size();
  int n_clusters = M_hat_.size();
  int n_gamma = pi_.size();
  NumericMatrix u_star_(n_clusters, n_gamma);
  
  for (int i = 0; i < n_clusters; ++i) {
    NumericMatrix M_hat_i = M_hat_(i);
    
    for (int r = 0; r < n_gamma; ++r) {
      NumericVector pi_r = pi_(r);
      
      double integration = 0;
      for (int s = 1; s < n_timesteps; ++s) {
        
        double inner_sum = 0;
        for (int j = 0; j < M_hat_i.nrow(); ++j) {
          inner_sum += M_hat_i(j, s) - M_hat_i(j, s - 1);
        }
        integration += pi_r(s) * p_hat(s - 1) * inner_sum/Ycal_(s);
      }
      
      u_star_(i,r) = integration;
    }
  }
  
  return u_star_;
}
