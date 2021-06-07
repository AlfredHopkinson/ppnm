

#ifndef HAVE_MINI_H
#define HAVE_MINI_H


void gradient(double F(gsl_vector * x), gsl_vector * x, gsl_vector * grad);

double backtrack(void (*f)(gsl_vector *x, double * fx), gsl_vector * x, gsl_vector * gf, gsl_vector * dx, gsl_matrix *H, int type);

int quasinewton(double F(gsl_vector*x), gsl_vector *x, double eps);


#endif


