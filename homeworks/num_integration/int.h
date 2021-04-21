#ifndef HAVE_int_H
#define HAVE_int_H

void integrate(double (*f)(double), double a, double b, double d, double e, double f2, double f3, int nrec);
void recadapt(double (*f)(double), double a, double b, double d, double e);


#endif
