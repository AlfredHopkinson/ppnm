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
#include"oldfunc.h"


















//maked the exaple test x^4 here

void example(gsl_vector *x, gsl_vector * fx){
	double t, e;
	t = gsl_vector_get(x,0);
	e = t*t*t;
	gsl_vector_set(fx,0,e);
}





//make newton func from the basics in the chapter
//using the hjacobi so pull a lot from previous homeworks

void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double eps){
	double d_x = sqrt(DBL_EPSILON); //use mnachine epsilon
	int n = x->size;


	gsl_vector * fx = gsl_vector_alloc(n);
	gsl_vector * dx = gsl_vector_alloc(n);
	gsl_vector * df = gsl_vector_alloc(n);
	gsl_matrix * J = gsl_matrix_alloc(n,n); 
	gsl_matrix * R = gsl_matrix_alloc(n,n);
	gsl_vector * y = gsl_vector_alloc(n);
	gsl_vector * fy = gsl_vector_alloc(n);

	f(x,fx);
	for(int j=0; j<n; j++){
		gsl_vector_set(x,j,gsl_vector_get(x,j)+d_x);
		f(x,df);
		gsl_vector_sub(df,fx);

		for(int i=0;i<n;i++){
			gsl_matrix_set(J,i,j,gsl_vector_get(df,i)/d_x);
			gsl_vector_set(x,i,gsl_vector_get(x,i)-d_x);
		}
		gsl_vector_scale(fx,-1);
		GS_decomp(J,R);
		GS_solve(J,R,fx,dx);

		double s = 2;

		while (1){ //here the python code example wanted true so i have put 1 
			s /= 2;
			gsl_vector_memcpy(y,x);
			gsl_blas_daxpy(s,dx,y);
			f(y,fy);
			if ((gsl_blas_dnrm2(fy)<(1.0-s/2)*gsl_blas_dnrm2(fx)) || (s<1/64)) break;
		}
		gsl_vector_memcpy(fx,fy);
		if (gsl_blas_dnrm2(fx)<eps||gsl_blas_dnrm2(dx)<d_x) break;
	}
	//forgot to free all the vectors 
	gsl_vector_free(fx);
	gsl_vector_free(dx);
	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(y);
	gsl_vector_free(fy);
}

int main(){
	//start with simplke to debug
	printf("Simple check to start\n");
	printf("\n");
	//not to self the newton needs f, x, eps
	gsl_vector * x = gsl_vector_alloc(1);
	double eps = 0.00001;
	//we also need a value to be in the x vector
	gsl_vector_set(x,0,5);

	newton(example,x,eps);
	printf("the root of x^4 is %10g\n",gsl_vector_get(x,0));

return 0;
}





