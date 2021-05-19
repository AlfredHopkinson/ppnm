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


void rkstep12(
	void f(double t,gsl_vector*y,gsl_vector*dydt), /* the f from dy/dt=f(t,y) */
	double t,              /* the current value of the variable */
	gsl_vector* yt,            /* the current value y(t) of the sought function */
	double h,              /* the step to be taken */
	gsl_vector* yh,             /* output: y(t+h) */
	gsl_vector* err){
	int c = yt->size;
	gsl_vector* k0 = gsl_vector_alloc(c);
	gsl_vector* k12 = gsl_vector_alloc(c);
	gsl_vector* ytt = gsl_vector_alloc(c);
	
	f(t,yt,k0); // this takes the step and evaluates 
	for(int i = 0; i<c; i++){
		gsl_vector_set(ytt,i, gsl_vector_get(yt,i)+gsl_vector_get(k0,i)*0.5*h);
	}
	f(t+0.5*h,ytt,k12);
	for(int i=0; i<c; i++){
		gsl_vector_set(yh,i, gsl_vector_get(yt,i)+gsl_vector_get(k12,i)*h);
	}
	for(int i=0;i<c;i++){
		gsl_vector_set(err,i, (gsl_vector_get(k0,i)-gsl_vector_get(k12,i))*0.5);
	}
	gsl_vector_free(k0);
	gsl_vector_free(k12);
	gsl_vector_free(ytt);
}


int driver(
	void (*f)(double t,gsl_vector* y,gsl_vector* dydt), /* right-hand-side of dy/dt=f(t,y) */
	double a,                     /* the start-point a */
	gsl_vector* ya,                     /* y(a) */
	double b,                     /* the end-point of the integration */
	gsl_vector* yb,                     /* y(b) to be calculated */
	double h,                     /* initial step-size */
	double acc,                   /* absolute accuracy goal */
	double eps,
	gsl_vector* err
												
){
	int k = ya->size;
	double t = a; //this will be the starting point  coming back to this added in a parint to show the starting point want working take out later
	printf("%10g ",t);
	for(int i=0;i<k;i++){
		printf("%10g ",gsl_vector_get(ya,i));
	}
	printf("\n");
	while(t<b){
		if(t+h>b) h = b-t;//a check that it doesnt step too far
			rkstep12(f,t,ya,h,yb,err); //put it through the RK with the y(t=h) being now in yb

			double nerr = gsl_blas_dnrm2(err);
			double normyb = gsl_blas_dnrm2(yb);
			double tol = (eps*normyb+acc)*sqrt(h/(b-a));
			
			if(nerr<tol){
				t = t+h;
				gsl_vector_memcpy(ya,yb);
				printf("%10g ",t);
				for(int i=0;i<k;i++){
					printf("%10g ",gsl_vector_get(ya,i));
				}
				printf("\n");
			}
			h *= pow(tol/nerr,0.25)*0.95; 

	}
}


void GS_decomp(gsl_matrix *A, gsl_matrix *R){
	int m = A->size2;
	int n = A->size1;
	assert(A->size2 == R ->size1);
	double Rij, Rmatrixii;
	for(int i=0;i<m;i++){
		gsl_vector_view icolumn = gsl_matrix_column(A, i);
		Rmatrixii= gsl_blas_dnrm2(&icolumn.vector); 
		gsl_matrix_set(R, i, i, Rmatrixii);
		gsl_vector_scale(&icolumn.vector,1.0/Rmatrixii);

		for(int j=i+1;j<m;j++){
			gsl_vector_view columnj = gsl_matrix_column(A,j);
			gsl_blas_ddot(&icolumn.vector, &columnj.vector, &Rij);	
			gsl_matrix_set(R,i,j,Rij); 
			gsl_blas_daxpy(-Rij, &icolumn.vector, &columnj.vector);
		}
	}
}


void GS_solve(gsl_matrix* Qb, gsl_matrix* Rb, gsl_vector* b, gsl_vector* x){
	gsl_blas_dgemv(CblasTrans,1,Qb,b,0,x);
	for(int i=x->size-1;i>=0;i--){
		double s=gsl_vector_get(x,i);
		for(int k=i+1;k<x->size;k++){
			s-=gsl_matrix_get(Rb,i,k)*gsl_vector_get(x,k);
		}
		gsl_vector_set(x,i,s/gsl_matrix_get(Rb,i,i));
						
	}
}










