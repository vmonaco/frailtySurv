#include <Rcpp.h>

typedef double (*density_fn)(double, double*);
typedef double (*deriv_density_fn)(double, double*, int);
typedef double (*deriv_deriv_density_fn)(double, double*, int, int);
typedef double (*lt_fn)(int, double, double*);

double lt_dgamma(int, double, double*);
double lt_dpvf(int, double, double*);
double lt_dlognormal(int, double, double*, double, double, int);
double lt_dinvgauss(int, double, double*, double, double, int);

double deriv_lt_dgamma(int, double, double*, int);
double deriv_lt_dpvf(int, double, double*, int);
double deriv_lt_dlognormal(int, double, double*, int, double, double, int);
double deriv_lt_dinvgauss(int, double, double*, int, double, double, int);

double deriv_deriv_lt_dgamma(int, double, double*, int, int);
double deriv_deriv_lt_dpvf(int, double, double*, int, int);
double deriv_deriv_lt_dlognormal(int, double, double*, int, int, double, double, int);
double deriv_deriv_lt_dinvgauss(int, double, double*, int, int, double, double, int);
