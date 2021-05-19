#ifndef HAVE_oldfunc_H
#define HAVE_oldfunc_H

void print_vector(const char* p, gsl_vector * A);
void print_matrix(const char* p, gsl_matrix * A);
void rkstep12(void f(double t,gsl_vector*y,gsl_vector*dydt), double t, gsl_vector* yt, double h, gsl_vector* yh, gsl_vector* err);
void driver(void (*f)(double t,gsl_vector* y,gsl_vector* dydt),double a, gsl_vector* ya, double b, gsl_vector* yb,  double h, double acc,  double eps,gsl_vector* err);



#endif
