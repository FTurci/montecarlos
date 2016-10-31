#ifndef __COMPASS_SEARCH__
#define __COMPASS_SEARCH__

#include <string>
#include "liquid.h"
double *compass_search ( double function_handle ( int m, double x[] ), int m, 
  double x0[], double delta_tol, double delta_init, int k_max, double &fx, 
  int &k );
double r8_abs ( double x );
void r8vec_copy ( int n, double a1[], double a2[] );
void r8vec_print ( int n, double a[], std::string title );
void timestamp ( );
 double  my_f ( int N3,double X[], liquid M);

#endif
