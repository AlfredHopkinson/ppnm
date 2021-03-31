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
	
	double n = 4;
	//create the random matrix A from my previous code
	gsl_matrix* A=gsl_matrix_alloc(n,n);
	gsl_matrix* D=gsl_matrix_alloc(n,n);
	gsl_matrix* V=gsl_matrix_alloc(n,n);
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++){
			gsl_matrix_set(A,i,j,rand()%10);
			gsl_matrix_set(V,i,j,0); //had a bit of trouble here as needed to make V a diagonal
		}
		gsl_matrix_set(V,i,i,1);
	}
	print_matrix("A=",A);
	//make a copy
	gsl_matrix_memcpy(D,A);
	//now I have A and V i can do the diagonalisation 
	print_matrix("V=",V);
	jacobi_diag(D,V);
	print_matrix("D=",D);
	print_matrix("V=",V);


return 0;
}

