#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#define RND (double)rand()/RAND_MAX
#ifndef NDEBUG
#define TRACE(args...) fprintf(stderr,args)
#else
#define TRACE(...)
#endif

void vector_print(char s[], gsl_vector* v){
	printf("%s\n",s);
	for(int i=0;i<v->size;i++)printf("%10g ",gsl_vector_get(v,i));
	printf("\n");
}



int main(){
	int n=5;
TRACE("n=%i\n",n);
TRACE("1\n");
	gsl_matrix* A=gsl_matrix_alloc(n,n);
	gsl_matrix* Acopy=gsl_matrix_alloc(n,n);
	gsl_vector* b=gsl_vector_alloc(n);
	gsl_vector* x=gsl_vector_alloc(n);
	gsl_vector* y=gsl_vector_calloc(n);
TRACE("2\n");
	for(int i=0; i< A->size1; i++)
		for(int j=0; j<A->size2; j++)
		{
		double Aij=RND;
		gsl_matrix_set(A,i,j,Aij);
		}
TRACE("3\n");
	gsl_matrix_memcpy(Acopy,A);
TRACE("4\n");
	for(int i=0; i< b->size; i++)
		{
		double bi=RND;
		gsl_vector_set(b,i,bi);
		}
TRACE("5\n");
	gsl_linalg_HH_solve(Acopy,b,x);
TRACE("6\n");
	gsl_blas_dgemv(CblasNoTrans,1,A,x,0,y);
TRACE("7\n");
	vector_print("right-hand side b:",b);
	vector_print("check: A*x should be equal b:",y);

TRACE("8\n");
gsl_matrix_free(A);
gsl_matrix_free(Acopy);
gsl_vector_free(b);
gsl_vector_free(x);
gsl_vector_free(y);
TRACE("9\n");
return 0;
}



//can use gdb main to debug but easier to printf debug
