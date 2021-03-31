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
			printf("%10g",gsl_matrix_get(A,i,j));
		}
		printf("\n");
	}
}

//something like in the example 

void timesJ(gsl_matrix* A, int p, int q, double omega){
	double c=cos(omega),s=sin(omega);
	for(int i=0;i<A->size1;i++){
		double newA_ip=c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
		double newA_iq=s*gsl_matrix_get(A,i,p)+c*gsl_matrix_get(A,i,q);
		gsl_matrix_set(A,i,p,newA_ip);
		gsl_matrix_set(A,i,q,newA_iq);
	}
}

void Jtimes(gsl_matrix* A, int p, int q, double omega){
	double c=cos(omega),s=sin(omega);
	for(int j=0;j<A->size2;j++){
		double newA_pj= c*gsl_matrix_get(A,p,j)+s*gsl_matrix_get(A,q,j);
		double newA_qj=-s*gsl_matrix_get(A,p,j)+c*gsl_matrix_get(A,q,j);
		gsl_matrix_set(A,p,j,newA_pj);
		gsl_matrix_set(A,q,j,newA_qj);
	}
}

void jacobi_diag (gsl_matrix * A, gsl_matrix * V){
	int changed;
	int n = A->size1;
	double omega , c, s;
	double apq, aqq, app, new_app, new_aqq;

	do{
		changed=0;
		for(int p=0;p<n-1;p++){
			for(int q=p+1;q<n;q++){

				apq=gsl_matrix_get(A,p,q);
				app=gsl_matrix_get(A,p,p);
				aqq=gsl_matrix_get(A,q,q);
				omega=0.5*atan(2*apq/(aqq-app));
				c=cos(omega),s=sin(omega);
		
				new_app=c*c*app-2*s*c*apq+s*s*aqq;
				new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
				if(new_app!=app || new_aqq!=aqq) // do rotation
					{
					changed=1;
					timesJ(A,p,q, omega);
					Jtimes(A,p,q,-omega); // A←J^T*A*J 
					timesJ(V,p,q, omega); // V←V*J
				}
			}
		}
	}while(changed!=0);
}







