#ifndef HAVE_leastsf_H
#define HAVE_leastsf_H

void print_vector(const char* p, gsl_vector * A);
void print_matrix(const char* p, gsl_matrix * A);
void rkstep12(void f(double t,gsl_vector*y,gsl_vector*dydt) /* the f from dy/dt=f(t,y) */
	double t,              /* the current value of the variable */
	gsl_vector* yt,            /* the current value y(t) of the sought function */
	double h,              /* the step to be taken */
	gsl_vector* yh,             /* output: y(t+h) */
	gsl_vector* err)
void driver(void(double t,vector* y,vector* dydt), /* right-hand-side of dy/dt=f(t,y) */
	double a,                     /* the start-point a */
	gsl_vector* ya,                     /* y(a) */
	double b,                     /* the end-point of the integration */
	gsl_vector* yb,                     /* y(b) to be calculated */
	double h,                     /* initial step-size */
	double acc,                   /* absolute accuracy goal */
	double eps,
	double s,
	double s1,
	double err,
	double normyb,
	double tol)
#endif
