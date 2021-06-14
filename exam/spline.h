#include"gsl/gsl_vector.h"
#ifndef HAVE_ANN_H
#define HAVE_ANN_H

int binsearch(int n, double* x, double z);
typedef struct {int n; double *x,*y,*c,*d;} akima_spline;
akima_spline* akima_spline_alloc(int n, double *x, double *y);
double akima_spline_eval(akima_spline *s, double z);
void akima_spline_free(akima_spline *s);

#endif
