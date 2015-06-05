#include <Rcpp.h>

typedef double (*density_fn)(double, double*);
typedef double (*deriv_density_fn)(double, double*, int);

struct phi_params {
  int k;
  int N_dot;
  double H_dot;
  double *theta;
  density_fn density;
};

struct phi_prime_params {
  int k;
  int N_dot;
  double H_dot;
  double *theta;
  deriv_density_fn deriv_density;
  int deriv_idx;
};

double dgamma(double, double*);
double lt_dgamma(int, double, double*);
double deriv_dgamma(double, double*, int);
double phi(double, void*);