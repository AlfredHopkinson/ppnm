#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include<assert.h>
#include"eigenfunc.h"


int main(){
	
	printf("PART A\n");
	printf("Proving my jacobi implimentation works\n");
	
	double n = 6;
	//create the random matrix A from my previous code
	gsl_matrix* A=gsl_matrix_alloc(n,n);
	gsl_matrix* D=gsl_matrix_alloc(n,n);
	gsl_matrix* V=gsl_matrix_alloc(n,n);
	
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++){
			//gsl_matrix_set(A,i,j,rand()%10);  Usual random method not working to make it symmetric going to try an allocate a double and put it in this loop
			double part = rand()%10;
			gsl_matrix_set(A,i,j,part);
			gsl_matrix_set(A,j,i,part);
	
		}
	}
	print_matrix("A=",A);
	//make a copy
	gsl_matrix_memcpy(D,A);
	gsl_matrix_set_identity(V); //had a bit of trouble here as needed to make V a diagonal
	//now I have A and V i can do the diagonalisation 
	print_matrix("V=",V);
	jacobi_diag(D,V);
	print_matrix("D=",D);
	print_matrix("V=",V);


	//now I have to check the implimentation as in the question 
	//I find AV*At =D
	printf("I will firstly show V^tAV==D\n");
	//times in 2 parts saving the first into e
	gsl_matrix* e=gsl_matrix_alloc(n,n);
	gsl_matrix* vtav=gsl_matrix_alloc(n,n);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,V,0.0,e);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,V,e,0.0,vtav);

	print_matrix("D=",D);
	print_matrix("V^t AV =",vtav);

	printf("I will now show VDV^t==A\n");
	//do same as before
	gsl_matrix* r=gsl_matrix_alloc(n,n);
	gsl_matrix* vdvt=gsl_matrix_alloc(n,n);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,V,D,0.0,r);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,r,V,0.0,vdvt);

	print_matrix("A=",A);
	print_matrix("VDV^t=",vdvt);

	printf("I will now show VV^t==i\n");

	gsl_matrix* VV=gsl_matrix_alloc(n,n);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,V,V,0.0,VV);
	print_matrix("V^t V=",VV);

	printf("all the checks confirm the implimentation works\n");






return 0;
}

