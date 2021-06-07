

#ifndef HAVE_MINI_H
#define HAVE_MINI_H

void gradient(void (*f)(gsl_vector * x, double * df), gsl_vector * x, gsl_vector * grad);

double backtrack(double f(gsl_vector *x), gsl_vector * x, gsl_vector * gf, gsl_vector * dx, gsl_matrix *H, int type);

int quasinewton(double f(gsl_vector *x),void gradient(gsl_vector *x, gsl_vector *df),gsl_vector *x,double eps);


#endif


