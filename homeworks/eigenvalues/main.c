#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include<assert.h>
#include"eigenfunc.h"


//void print_matrix(const char* p, gsl_matrix * A){
//	printf("%s\n",p);
//	for(int i=0;i<A->size1;i++){
//		for(int j=0;j<A->size2;j++){
//			printf("%10g",gsl_matrix_get(A,i,j));
//		}
//		printf("\n");
//	}
//}


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

	//
	//Part B starts here
	
	printf("\n\n\n");
	printf("Part B starts here\n\n");

	//build the hamiltonian matrix
	
	int nn = 100;//started with 20 but this made an awful graph so I will crank it up
	double s=1.0/(nn+1);
	gsl_matrix* H = gsl_matrix_alloc(nn,nn);
	for(int i = 0;i<nn-1;i++){
		gsl_matrix_set(H,i,i,-2);
		gsl_matrix_set(H,i,i+1,1);
		gsl_matrix_set(H,i+1,i,1);
	}
	gsl_matrix_set(H,nn-1,nn-1,-2);
	gsl_matrix_scale(H,-1/s/s);

	//diagonalize with my jacobhi
	
	gsl_matrix* G = gsl_matrix_alloc(nn,nn);
	jacobi_diag(H,G);

//	print_matrix("G=",G);



	//check energies are correct
	FILE * energies = fopen("energyout.txt","w");	


	for (int k=0;k<nn/3;k++){
		double exact = M_PI*M_PI*(k+1)*(k+1);
		double calculated = gsl_matrix_get(H,k,k);
		fprintf(energies,"%i %g %g\n", k, calculated, exact);
	}
	fclose(energies);

	int flip[3] = {-1,-1,-1};
	double f = 1.0/sqrt(s);



	FILE * quantuneigenvalues = fopen("quantumeigenvalues.txt","w");
	//this is all well and good but I cant get it to plot only one set at a time so I am going to alter it so I can
	fprintf(quantuneigenvalues, "%g %g %g %g\n",0.0,0.0,0.0,0.0);
	for (int k=0; k<nn;k++){
//		fprintf(quantuneigenvalues,"%g %g\n",0.0,0.0);
//		for(int i=0; i<nn; i++)
//			fprintf(quantuneigenvalues,"%g %g\n",(i+1.0)/(nn+1), gsl_matrix_get(G,i,k));
//		fprintf(quantuneigenvalues,"%g %g\n",1.0,0.0);
		fprintf(quantuneigenvalues, "%g %g %g %g\n",(k+1.0)/(nn+1),-1.0*gsl_matrix_get(G,k,0),gsl_matrix_get(G,k,1),gsl_matrix_get(G,k,2));
	}
	fprintf(quantuneigenvalues, "%g %g %g %g\n",1.0,0.0,0.0,0.0);
	fclose(quantuneigenvalues);

	gsl_matrix_free(H);
	gsl_matrix_free(G);





return 0;
}

