#ifndef oldfunc_H
#define oldfunc_H

void print_vector(const char* p, gsl_vector * A);
void print_matrix(const char* p, gsl_matrix * A);
void rkstep12(void f(double t,gsl_vector*y,gsl_vector*dydt), double t, gsl_vector* yt, double h, gsl_vector* yh, gsl_vector* err);
void driver(void (*f)(double t,gsl_vector* y,gsl_vector* dydt),double a, gsl_vector* ya, double b, gsl_vector* yb,  double h, double acc,  double eps,gsl_vector* err);
void GS_decomp(gsl_matrix *A, gsl_matrix *R);
void GS_solve(gsl_matrix* Qb, gsl_matrix* Rb, gsl_vector* b, gsl_vector* x);


#endif
