#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include<assert.h>
#include <gsl/gsl_integration.h>
#define RND ((double)rand()/RAND_MAX)
#include<float.h>




//we need to make the quasinewton func and have the implimentatiopn template to work off.
double DELTA = sqrt(DBL_EPSILON);
void qnewton(void (*f)(gsl_vector * x, gsl_vector * fx),gsl_vector * x, double eps){
	int n = x->size, numsteps = 0;
	
	gsl_vector * gfx = gsl_vector_alloc(n);
	gsl_vector * s = gsl_vector_alloc(n);

	while (numsteps<1000){
		if(gsl_blas_dnrm2(gfx)<eps*eps) break;
		gsl_blas_dgemv(CblasNoTrans,-1,B,gfx,0,e);


















