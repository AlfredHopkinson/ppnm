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





