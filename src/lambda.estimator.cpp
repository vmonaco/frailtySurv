#include <RcppGSL.h>
// [[Rcpp::depends(RcppGSL)]]

#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_psi.h>

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

// [[Rcpp::export]]
NumericVector dgamma_(NumericVector x, double theta) {
  int n = x.size();
  NumericVector out(n);
  for(int i = 0; i < n; ++i) {
    out[i] = gsl_ran_gamma_pdf(x[i], 1/theta, theta);
  }
  return out;
}

// [[Rcpp::export]]
double deriv_dgamma(double x, double theta) {
//   if (params.size() != 1) {
//     throw std::range_error("params must be length 1 for the gamma distribution");
//   }
//   
//   // deriv_idx is an R index for the parameter to take a partial derivative against
//   if (deriv_idx <= 0 || deriv_idx > params.size()) {
//     throw std::range_error("deriv_idx out of range, must be > 0 and <= length(params)");
//   }
  
//   int n = x.size();
//   NumericVector out(n);
  double term1, term2, term3;
  // for(int i = 0; i < n; ++i) {
      term1 = pow((x/theta),(1/theta - 1));
      term2 = exp(-x/theta);
      term3 = log(theta/x) + gsl_sf_psi(1/theta) + x - 1;
      double out = (term1 * term2 * term3)/(tgamma(1/theta)*pow(theta,3));
  // }
  return out;
}

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

double phi_1(double w, void* data) {
  double (*params)[3] = (double (*)[3]) data;
  double N_dot = (*params)[0];
  double H_dot = (*params)[1];
  double theta = (*params)[2];
  return pow(w, N_dot) * exp(-w * H_dot) * gsl_ran_gamma_pdf(w, 1/theta, theta);
}

double phi_1_deriv(double w, void* data) {
  double (*params)[3] = (double (*)[3]) data;
  double N_dot = (*params)[0];
  double H_dot = (*params)[1];
  double theta = (*params)[2];
  return pow(w, N_dot) * exp(-w * H_dot) * deriv_dgamma(w, theta);
}

double phi_2(double w, void* data) {
  double (*params)[3] = (double (*)[3]) data;
  double N_dot = (*params)[0];
  double H_dot = (*params)[1];
  double theta = (*params)[2];
  return pow(w, N_dot + 1) * exp(-w * H_dot) * gsl_ran_gamma_pdf(w, 1/theta, theta);
}

// [[Rcpp::export]]
double psi(double theta, int N_dot, double H_dot, String frailty_distr) {
  
  if (frailty_distr == "gamma") {
    // gamma shortcut
    // return (N_dot + 1/theta)/(H_dot + 1/theta);
  } else {
    throw std::range_error("Unsupported frailty distribution");
  }
  
  double data[3];
  data[0] = N_dot;
  data[1] = H_dot;
  data[2] = theta;
  double out = integrate(&phi_2, &data)/integrate(&phi_1, &data);
  return out;
}

// [[Rcpp::export]]
Rcpp::List baseline_hazard_estimator(Rcpp::List X_, Rcpp::List k_,
                                     NumericVector d_,
                                     Rcpp::List Y_, Rcpp::List N_dot,
                                     NumericVector beta, double theta, 
                                     String frailty_distr) {
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
      
      denom += tmp * psi(theta, N_dot_i[k-1], H_dot_i[k-1], frailty_distr);
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
double U_r(Rcpp::List R_, Rcpp::List N_dot,
               Rcpp::List I_, Rcpp::List X_,
               Rcpp::List H_, Rcpp::List H_dot,
               NumericVector beta, double theta, 
               int beta_idx, String frailty_distr) {
  
  int n_clusters = X_.size();
  
  double term1 = 0;
  double term2 = 0;
  for (int i = 0; i < n_clusters; ++i) {
    Rcpp::NumericMatrix X_i = X_[i];
    Rcpp::NumericMatrix H_i = H_[i];
    Rcpp::NumericVector I_i = I_[i]; // Failure indicator
    Rcpp::NumericVector R_i = R_[i]; // Failure rank, whare T[r] = failure time
    Rcpp::NumericVector N_dot_i = N_dot[i];
    Rcpp::NumericVector H_dot_i = H_dot[i];
    
    double tmp = 0;
    for (int j = 0; j < X_i.nrow(); j++) {
      term1 += I_i[j] * X_i(j,beta_idx - 1);
      tmp += H_i(j, R_i[j] - 1) * X_i(j,beta_idx - 1);
//       Rcout << "term1: " << term1 << std::endl;
//       Rcout << "tmp: " << tmp << std::endl;
//       Rcout << "H_i: " << H_i(j, R_i[j] - 1) << std::endl;
//       Rcout << "R_i: " << R_i[j] << std::endl;
//       Rcout << "X_i: " << X_i(j, beta_idx-1) << std::endl;
      
    }
//     Rcout << "N_i: " << N_dot_i[N_dot_i.size()-1] << std::endl;
//     Rcout << "psi: " << psi(theta, N_dot_i[N_dot_i.size()-1], H_dot_i[H_dot_i.size()-1], "gamma") << std::endl;
    term2 += tmp * psi(theta, N_dot_i[N_dot_i.size()-1], H_dot_i[H_dot_i.size()-1], "gamma");
  }
//   Rcout << "term1: " << term1 << std::endl;
//   Rcout << "term2: " << term2 << std::endl;
  term1 = term1/n_clusters;
  term2 = term2/n_clusters;
  return term1 - term2;
}

// [[Rcpp::export]]
double U_p(Rcpp::List R_, Rcpp::List N_dot,
           Rcpp::List I_, Rcpp::List X_,
           Rcpp::List H_, Rcpp::List H_dot,
           NumericVector beta, double theta,
           int theta_idx, int tau_k, String frailty_distr) {
  
  int n_clusters = X_.size();
  
  double out = 0;
  for (int i = 0; i < n_clusters; ++i) {
    Rcpp::NumericMatrix X_i = X_[i];
    Rcpp::NumericMatrix H_i = H_[i];
    Rcpp::NumericVector I_i = I_[i]; // Failure indicator
    Rcpp::NumericVector R_i = R_[i];
    Rcpp::NumericVector N_dot_i = N_dot[i];
    Rcpp::NumericVector H_dot_i = H_dot[i];
    
    double data[3];
    data[0] = N_dot_i[tau_k - 1];
    data[1] = H_dot_i[tau_k - 1];
    data[2] = theta;
    double numer = integrate(&phi_1_deriv, &data);
    double denom = integrate(&phi_1, &data);
    out += numer/denom;
  }
  
  out = out/n_clusters;
  
  return out;
}

// [[Rcpp::export]]
double log_likelihood(NumericVector beta, double theta, 
                      NumericVector lambda, List H_dot,
                      List I_, List R_, List X_, List N_dot,
                      String frailty_distr) {
  
  int n_clusters = X_.size();
  
  double term1 = 0;
  double term2 = 0;
  for (int i = 0; i < n_clusters; ++i) {
    Rcpp::NumericMatrix X_i = X_[i];
    Rcpp::NumericVector I_i = I_[i]; // Failure indicator
    Rcpp::NumericVector R_i = R_[i]; // Failure rank, whare T[r] = failure time
    Rcpp::NumericVector N_dot_i = N_dot[i];
    Rcpp::NumericVector H_dot_i = H_dot[i];
    
    double term1_tmp = 0;
    for (int j = 0; j < X_i.nrow(); j++) {
      term1_tmp += I_i(j) * log(lambda(R_i[j] - 1) * exp(sum(beta * X_i[j])));
    }
    
    term1 += term1_tmp;
    
    double data[3];
    data[0] = N_dot_i[N_dot_i.size() - 1];
    data[1] = H_dot_i[H_dot_i.size() - 1];
    data[2] = theta;
    term2 += log(integrate(&phi_1, &data));
  }
  
  return term1 + term2;
}