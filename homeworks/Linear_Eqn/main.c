#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>










int main(){

	// Part A Convert to QR and then apply gram-schmidt
	printf("PART A\n");
	printf("Creating the function and performing the Gram-Schmidt on it\n");
	//size of the matrix  REMEMBER  it isnt square so n>m
	double n =6, m=5;
	//create the starting matrix and generate a random matrix 
	gsl_matrix* A=gsl_matrix_alloc(n,m);
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			gsl_matrix_set(A,i,j,100.0*rand()/RAND_MAX);
		}
	}

	//Factorisation need to make both Q and R and use Gram
	
	gsl_matrix* Q=gsl_matrix_alloc(n,m);
	gsl_matrix* R=gsl_matrix_alloc(m,m);

	GS_decomp(Q,R);





return 0;
}
