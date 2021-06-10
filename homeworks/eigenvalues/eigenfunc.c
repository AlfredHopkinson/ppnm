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

//something like in the example 

void timesJ(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int i=0;i<A->size1;i++){
		double new_aip=c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
		double new_aiq=s*gsl_matrix_get(A,i,p)+c*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,new_aip);
		gsl_matrix_set(A,i,q,new_aiq);
	}
}

void Jtimes(gsl_matrix* A, int p, int q, double theta){
	double c=cos(theta),s=sin(theta);
	for(int j=0;j<A->size1;j++){
		double new_apj= c*gsl_matrix_get(A,p,j)+s*gsl_matrix_get(A,q,j);
		double new_aqj=-s*gsl_matrix_get(A,p,j)+c*gsl_matrix_get(A,q,j);
		gsl_matrix_set(A,p,j,new_apj);
		gsl_matrix_set(A,q,j,new_aqj);
	}
}

void jacobi_diag (gsl_matrix * A, gsl_matrix * V){
	assert(A->size1 == A->size2);
	int changed;
	int n = A->size1;
	gsl_matrix_set_identity(V);	
	
	do{
		changed=0;
		for(int p=0;p<A->size1;p++){
			for(int q=p+1;q<n;q++){
				double apq=gsl_matrix_get(A,p,q);
				double app=gsl_matrix_get(A,p,p);
				double aqq=gsl_matrix_get(A,q,q);
				double theta=0.5*atan2(2*apq,aqq-app);
				double c=cos(theta),s=sin(theta);
				double new_app=c*c*app-2*s*c*apq+s*s*aqq;
				double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
				if(new_app!=app || new_aqq!=aqq) 
					{
					changed=1;
					timesJ(A,p,q, theta);
					Jtimes(A,p,q,-theta); 
					timesJ(V,p,q, theta);
					}
				}
			}
	}while(changed!=0);
}







