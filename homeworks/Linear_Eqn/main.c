#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include<assert.h>

void GS_decomp(gsl_matrix *A, gsl_matrix *R){
	int m = A->size2;
	int n = A->size1;
	assert(A->size2 == R ->size1);
	assert(n > m);
	double Rij, Rmatrixii;
	for(int i=0;i<m;i++){
		gsl_vector_view icolumn = gsl_matrix_column(A, i);
		
		// the linear says we need to generate the equation but need to get the sqrt of a.a 
		// coulnt find a nice way to manually do this so found a blas
		Rmatrixii= gsl_blas_dnrm2(&icolumn.vector); //this blas does the exact right thing we need it for
		gsl_matrix_set(R, i, i, Rmatrixii);
		//Once we have set the new matrix parts we normalize it
		gsl_vector_scale(&icolumn.vector,1.0/Rmatrixii);

		for(int j=i+1;j<m;j++){
			gsl_vector_view columnj = gsl_matrix_column(A,j);
			gsl_blas_ddot(&icolumn.vector, &columnj.vector, &Rij);	// computes the scalar product 
			gsl_matrix_set(R,i,j,Rij); // set the numbers into R
			//still have to orthoganalise this 
			gsl_blas_daxpy(-Rij, &icolumn.vector, &columnj.vector);//blas for computing y = cx+y linear equations but had to change previous stuff to vectors now altered the rest but cant check if it is done so may need to get a way to print the matrix
		}
	}
}

//print a matrix 
void print_matrix(const char* p, gsl_matrix * A){
	printf("%s\n",p);
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++){
			printf("%g",gsl_matrix_get(A,i,j));
		}
		printf("\n");
	}
	printf("\n");
}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
	//apply qt to vector b and save in x using blas again
	gsl_blas_dgemv(CblasTrans,1,Q,b,0,x);
	//preform the back subs from the woorkbooklet
	
	for(int i=x->size-1;i>=0;i--){
		double s=gsl_vector_get(x,i);
		for(int k=i+1;k<x->size;k++){
			s-=gsl_matrix_get(R,i,k)*gsl_vector_get(x,k);
			gsl_vector_set(x,i,s/gsl_matrix_get(R,i,i));
		}
	}
}







int main(){

	// Part A Convert to QR and then apply gram-schmidt
	printf("PART A\n");
	printf("Creating the function and performing the Gram-Schmidt on it\n");
	//size of the matrix  REMEMBER  it isnt square so n>m
	double n =6, m=5;
	//create the starting matrix and generate a random matrix 
	gsl_matrix* A=gsl_matrix_alloc(n,m);
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++){
			gsl_matrix_set(A,i,j,rand()%10);
		}
	}

	//Factorisation need to make both Q and R and use Gram
	
	gsl_matrix* Q=gsl_matrix_alloc(n,m);
	gsl_matrix* R=gsl_matrix_alloc(m,m);
	gsl_matrix_memcpy(Q, A);

	GS_decomp(Q,R);
	//got the decomposition now do the checks - involves printing matrix 
	print_matrix("A=",A);
	print_matrix("R=",R);
	print_matrix("Q=",Q);

	//doing this shows the R being triangular and that it works next is to do the transpose 
//	gsl_matrix* Qt=gsl_matrix_alloc(m,n);
//	gsl_matrix* Qtcopy=gsl_matrix_alloc(m,n);
	
//	gsl_matrix_transpose_memcpy(Qt,Q);
//	gsl_matrix_memcpy(Qtcopy,Qt);
//	gsl_matrix_mul_elements(Qt,Q);
//	print_matrix("Qt*Q=",Qt);
//      didnt work due to different sizes of martic and i can find another way so trying blas


	gsl_matrix* QtimesQt = gsl_matrix_alloc(m,m);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans, 1.0, Q,Q,0.0, QtimesQt );//this does it for you which would have been much easier from the start
	printf("We can show Qt*Q=1");
	print_matrix("Qt*Q=", QtimesQt);
	printf("I believe this is close to 0 and only isnt due to rounding errors\n");
	
	
	gsl_matrix * QR=gsl_matrix_alloc(n,m);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, Q,R,0.0, QR );
	printf("Checking QR=A\n");
	print_matrix("QR=",QR);


	printf("Now time for part 2 GS_solve\n");

	//put in the solve funct
	int nm=5;

	gsl_matrix* Ab=gsl_matrix_alloc(nm,nm);
	gsl_matrix* b=gsl_matrix_alloc(nm,nm);
	
	for(int i=0;i<Ab->size1;i++){
		for(int j=0;j<Ab->size2;j++){
			gsl_matrix_set(Ab,i,j,rand()%10);
		}
	}
	for(int i=0;i<b->size1;i++){
		for(int j=0;j<b->size2;j++){
			gsl_matrix_set(b,i,j,rand()%10);
		}
	}








return 0;
}
