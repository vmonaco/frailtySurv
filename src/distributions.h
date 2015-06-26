#include <Rcpp.h>

typedef double (*density_fn)(double, double*);
typedef double (*deriv_density_fn)(double, double*, int);
typedef double (*deriv_deriv_density_fn)(double, double*, int, int);
typedef double (*lt_fn)(int, double, double*);

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

struct phi_prime_prime_params {
  int k;
  int N_dot;
  double H_dot;
  double *theta;
  deriv_deriv_density_fn deriv_deriv_density;
  int deriv_idx_1;
  int deriv_idx_2;
};

struct lt_params {
  int p;
  double s;
  double *theta;
  lt_fn lt;
  int deriv_idx;
};

double dgamma(double, double*);
double deriv_dgamma(double, double*, int);

double lt_dgamma(int, double, double*);
double lt_dposstab(int, double, double*);
double lt_dpvf(int, double, double*);

double deriv_lt_dgamma(int, double, double*, int);
double deriv_lt_dposstab(int, double, double*, int);
double deriv_lt_dpvf(int, double, double*, int);

double deriv_deriv_lt_dgamma(int, double, double*, int, int);

double dlognormal(double, double*);
double deriv_dlognormal(double, double*, int);

double dinvgauss(double, double*);
double deriv_dinvgauss(double, double*, int);