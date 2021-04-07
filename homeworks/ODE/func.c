#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include<assert.h>

void print_vector(const char* p, gsl_vector * A){
	printf("%s\n",p);
	for(int i=0;i<A->size;i++){
		printf("%10g",gsl_vector_get(A,i));
	}
	printf("\n");
												
}
void print_matrix(const char* p, gsl_matrix * A){
	printf("%s\n",p);
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++){
			printf("%15g",gsl_matrix_get(A,i,j));
		}
		printf("\n");
	}
}

//make the RK stepper we go for 12 like in the example
void rkstep12(
	//compared to the exampel x=t, yx=yt, h=h, yh=yh, dy= err
	void f(double t,gsl_vector*y,gsl_vector*dydt) /* the f from dy/dt=f(t,y) */
	double t,              /* the current value of the variable */
	gsl_vector* yt,            /* the current value y(t) of the sought function */
	double h,              /* the step to be taken */
	gsl_vector* yh,             /* output: y(t+h) */
	gsl_vector* err){

	int c = yt->size;
	//dont need to alloc yt as its called above and alloced in main
	gsl_vector* k0 = gsl_vector_alloc(c);
	gsl_vector* k12 = gsl_vector_alloc(c);
	gsl_vector* ytt = gsl_vector_alloc(c);


	f(c,t,yt,k0); // this takes the step and evaluates 
	for(int i = 0; i<c; i++){
		gsl_vector_set(yt,i, gsl_vector_get(yt)+gsl_vector_get(k0,i)*0.5*h);
	}
	f(c,t+0.5*h,yh,k12);
	for(int i=0; i<c; i++){
		gsl_vector_set(yh,i, gsl_vector_get(yt)+gsl_vector_get(k12,i)*h);
	}
	for(int i=0;i<c;i++){
		gsl_vector_set(err,i, (gsl_vector_get(k0,i)-gsl_vector_get(k12,i))*0.5*h;
	}
	gsl_vector_free(k0);
	gsl_vector_free(k12);

}

//now time to do the second driver part

void driver(
	void f(double,vector*,vector*), /* right-hand-side of dy/dt=f(t,y) */
	double a,                     /* the start-point a */
	vector* ya,                     /* y(a) */
	double b,                     /* the end-point of the integration */
	vector* yb,                     /* y(b) to be calculated */
	double h,                     /* initial step-size */
	double acc,                   /* absolute accuracy goal */
	double eps
){ 










