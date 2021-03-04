#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#define RND (double)rand()/RAND_MAX




void vector_print(char s[], gsl_vector* v){
 	printf("%s\n",s);
	for(int i=0;i<v->size;i++)printf("%10g ",gsl_vector_get(v,i));
	printf("\n");
}


int main(){
	int n=2;
	gsl_matrix* A=gsl_matrix_alloc(n,n);
	gsl_matrix* Acopy=gsl_matrix_alloc(n,n);
	gsl_vector* eval=gsl_vector_calloc(n);
	gsl_eigen_symm_workspace * W=gsl_eigen_symm_alloc(n); 	



	gsl_matrix_set(A,0,0,1.0);
	gsl_matrix_set(A,0,1,1.0/2);
	
	gsl_matrix_set(A,1,0,1.0/2);
	gsl_matrix_set(A,1,1,1.0/3);
	

	gsl_matrix_memcpy(Acopy,A);	
 																			
	gsl_eigen_symm(A,eval,W);


	vector_print("eigen values are =",eval);	


	gsl_matrix_free(A);
	gsl_matrix_free(Acopy);

return 0;
}
