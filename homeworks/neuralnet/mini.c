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
#include"mini.h"
#include"ann.h"

void gradient(void (*f)(gsl_vector * x, double * df), gsl_vector * x, gsl_vector * grad) {
	double fx, faux; int n = x->size;
	double DELTA = sqrt(DBL_EPSILON);
	f(x,&fx);
	for (int i = 0; i < n; i++) {
		gsl_vector_set(x,i,gsl_vector_get(x,i)+DELTA);
		f(x,&faux);
		faux -= fx;
		gsl_vector_set(grad,i,faux/DELTA);
		gsl_vector_set(x,i,gsl_vector_get(x,i)-DELTA);
	}
}


double backtrack(double f(gsl_vector *x), gsl_vector * x, gsl_vector * gf, gsl_vector * dx, gsl_matrix *H, int type){
	assert(type == 1 || type ==0);

	double l = 1, a = 0.0001, fx = f(x), deltaxdx;
	gsl_vector *xandl = gsl_vector_alloc(x->size);
	gsl_vector *ldx = gsl_vector_alloc((x->size));
					
	gsl_vector_memcpy(xandl,x);
	gsl_vector_memcpy(ldx,dx);
	gsl_vector_scale(ldx,l);
	gsl_vector_add(xandl,ldx);
	gsl_blas_ddot(dx,gf,&deltaxdx);
	while( f(xandl)>= fx +a * l * deltaxdx){
		l /= 2;
		gsl_vector_memcpy(xandl,x);
		gsl_vector_memcpy(ldx,dx);
		gsl_vector_scale(ldx,l);
		gsl_vector_add(xandl,ldx);

		if(type==1){
		if(l<0.000001){
		gsl_matrix_set_identity(H);
		break;}
		}
		else;
	}
	gsl_vector_free(xandl);
	gsl_vector_free(ldx);

	return l;
}



int quasinewton(double f(gsl_vector *x), gsl_vector *x,double eps){
	int n = x->size;
	int iter = 0;

	gsl_vector *df = gsl_vector_alloc(n);
	gsl_vector *y = gsl_vector_alloc(n);
	gsl_vector *s = gsl_vector_alloc(n);
	gsl_matrix *H = gsl_matrix_alloc(n,n);
	gsl_vector *Hy = gsl_vector_alloc(n);
	double HyTs;
	gsl_matrix_set_identity(H);

	gradient(x,df);

	while( gsl_blas_dnrm2(df)>eps){
		iter++;
		gsl_blas_dgemv(CblasNoTrans,-1.0,H,df,0.0,s);

		double lambda = backtrack(f,x,df,s,H,1);

		gsl_vector_scale(s,lambda);
		gsl_vector_add(x,s);
		gradient(x,y);
		gsl_vector_sub(y,df);
		gsl_blas_dgemv(CblasNoTrans,1.0,H,y,0.0,Hy);
		gsl_blas_dgemv(CblasNoTrans,1.0,H,s,0.0,df);
		gsl_blas_ddot(Hy,s,&HyTs);
		gsl_vector_sub(s,Hy);
		gsl_blas_dger(1.0/HyTs,s,df,H);

		gradient(x,df);
	}

	gsl_vector_free(df);
	gsl_vector_free(y);
	gsl_vector_free(s);
	gsl_matrix_free(H);
	gsl_vector_free(Hy);
	return iter;
}









