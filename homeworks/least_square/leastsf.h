#ifndef HAVE_leastsf_H
#define HAVE_leastsf_H

void print_vector(const char* p, gsl_vector * A);
void print_matrix(const char* p, gsl_matrix * A);
void GS_decomp(gsl_matrix *A, gsl_matrix *R);
void GS_solve(gsl_matrix* Qb, gsl_matrix* Rb, gsl_vector* b, gsl_vector* x);

#endif
