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



//maked the exaple test x^3 here

//void example(gsl_vector *x, gsl_vector * fx){
//	double t, e;
//	t = gsl_vector_get(x,0);
//	e = t*t*t;
//	gsl_vector_set(fx,0,e);
//}


void example(gsl_vector * x, gsl_vector * fx) {
	double x1 = gsl_vector_get(x,0);
	double x2 = gsl_vector_get(x,1);

	gsl_vector_set(fx,0,20*sin(x1)*cos(x2));
	gsl_vector_set(fx,1,20*cos(x2));
		
}


//Rosenbracks valley func here.
void rosenbrock (gsl_vector * R , gsl_vector * fR){
	double scale = 1000;
	double gx,gy;
	double x = gsl_vector_get(R,0);
	double y = gsl_vector_get(R,1);
	//now i worked out the gradients for both x and y 
	gx = -2*(1-x)-400*x*(y-(x*x));
	gy = 200*(y-(x*x));
	gsl_vector_set(fR,0,gx);
	gsl_vector_set(fR,1,gy);
	gsl_vector_scale(fR,scale);
}








//make newton func from the basics in the chapter
////using the hjacobi so pull a lot from previous homeworks

void newton(void f(gsl_vector* x,gsl_vector* fx), gsl_vector* x, double tol){
	double d_x = sqrt(DBL_EPSILON); //use mnachine epsilon
	int n = x->size;

	gsl_vector * fx = gsl_vector_alloc(n);
	gsl_vector * dx = gsl_vector_alloc(n);
	gsl_vector * df = gsl_vector_alloc(n);
	gsl_matrix * J = gsl_matrix_alloc(n,n); 
	gsl_matrix * R = gsl_matrix_alloc(n,n);
	gsl_vector * y = gsl_vector_alloc(n);
	gsl_vector * fy = gsl_vector_alloc(n);
	while (1){ //here the python code example wanted true so i have put 1
		f(x,fx);
		for(int j=0; j<n; j++){
			gsl_vector_set(x,j,gsl_vector_get(x,j)+d_x);
			f(x,df);
			gsl_vector_sub(df,fx);

			for(int i=0;i<n;i++){
				gsl_matrix_set(J,i,j,gsl_vector_get(df,i)/d_x);
				gsl_vector_set(x,i,gsl_vector_get(x,i)-d_x);
			}
		}
			gsl_vector_scale(fx,-1);
			GS_decomp(J,R);
			GS_solve(J,R,fx,dx);

			double s = 2;

			while (1){  
				s /= 2;
				gsl_vector_memcpy(y,x);
				gsl_blas_daxpy(s,dx,y);
				f(y,fy);
				if ((gsl_blas_dnrm2(fy)<(1.0-s/2)*gsl_blas_dnrm2(fx)) || (s<1/64)) break;
			}
			gsl_vector_memcpy(fx,fy);
			gsl_vector_memcpy(x,y);
			if (gsl_blas_dnrm2(fx)<tol||gsl_blas_dnrm2(dx)<d_x) break;
		}
	//forgot to free all the vectors 
	gsl_vector_free(fx);
	gsl_vector_free(dx);
	gsl_matrix_free(J);
	gsl_matrix_free(R);
	gsl_vector_free(y);
	gsl_vector_free(fy);
	gsl_vector_free(df);
}



int main(){
	//start with simplke to debug
	printf("Simple check to start\n");
	printf("\n");
	//not to self the newton needs f, x, eps
	gsl_vector * x = gsl_vector_alloc(2);
	double tol = 0.0001;
	//we also need a value to be in the x vector
	gsl_vector_set(x,0,10);

	newton(example,x,tol);
	printf("the root of x^4 is %10g\n\n",gsl_vector_get(x,0));
	gsl_vector_free(x);

	//now to do ropsenbrocks valley we need to make a func to contain and then pass it thought the newton again
	printf("The roots of Rosenbrock's valley function\n\n\n");
	tol = 0.001;
	gsl_vector * rosenvector = gsl_vector_alloc(2);
	double rosenx = 1.4;
	double roseny = 0.9;
//	gsl_vector_set(rosenvector,0,rosenx);
//	gsl_vector_set(rosenvector,1,roseny);
//	newton(rosenbrock,rosenvector,tol);
//	print_vector("The roots of the gradient of the Rosenbrock function equals = ",rosenvector);
	//wouldnt work but i allways forget to free the vectors afterwards
	gsl_vector_free(rosenvector);
	gsl_vector_free(x);







return 0;
}























