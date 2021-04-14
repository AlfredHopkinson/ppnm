#ifndef HAVE_leastsf_H
#define HAVE_leastsf_H

void print_vector(const char* p, gsl_vector * A);
void print_matrix(const char* p, gsl_matrix * A);
void rkstep12(void f(double t,gsl_vector*y,gsl_vector*dydt), /* the f from dy/dt=f(t,y) */
	double t,              /* the current value of the variable */
	gsl_vector* yt,            /* the current value y(t) of the sought function */
	double h,              /* the step to be taken */
	gsl_vector* yh,             /* output: y(t+h) */
	gsl_vector* err);

void driver(void f(double t, gsl_vector* y,gsl_vector* dydt),
	double a,
	gsl_vector* ya, 
	double b,
	gsl_vector* yb, 
	gsl_vector* err, 
	double h,
	double acc,
	double eps);






#endif
