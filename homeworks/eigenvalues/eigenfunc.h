#ifndef HAVE_leastsf_H
#define HAVE_leastsf_H

void print_vector(const char* p, gsl_vector * A);
void print_matrix(const char* p, gsl_matrix * A);
void timesJ(gsl_matrix* , int, int, double);
void Jtimes(gsl_matrix* , int, int, double);
void jacobi_diag(gsl_matrix*, gsl_matrix*);


#endif
