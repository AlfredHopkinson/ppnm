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
	//assert(n > m); had this but it breaks for part b so remember 
	double Rij, Rmatrixii;
	for(int i=0;i<m;i++){
		gsl_vector_view icolumn = gsl_matrix_column(A, i);
		
		// the tor("b =",b); says we need to generate the equation but need to get the sqrt of a.a 
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
			printf("%15g",gsl_matrix_get(A,i,j));
		}
		printf("\n");
	}
}

void GS_solve(gsl_matrix* Qb, gsl_matrix* Rb, gsl_vector* b, gsl_vector* x){
	//apply qt to vector b and save in x using blas again
	gsl_blas_dgemv(CblasTrans,1,Qb,b,0,x);
	//preform the back subs from the woorkbooklet
	
	for(int i=x->size-1;i>=0;i--){
		double s=gsl_vector_get(x,i);
		for(int k=i+1;k<x->size;k++){
			s-=gsl_matrix_get(Rb,i,k)*gsl_vector_get(x,k);
		}
		gsl_vector_set(x,i,s/gsl_matrix_get(Rb,i,i));
		
	}
}

void print_vector(const char* p, gsl_vector * A){
	printf("%s\n",p);
	for(int i=0;i<A->size;i++){
		printf("%10g",gsl_vector_get(A,i));
	}
	printf("\n");
	
}


//part b is to create an inverse function
//there is the basics of one on the chapter so modified that in here
void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
	//here there is the problem gs slove needs a vector x to put into
	gsl_vector*x = gsl_vector_alloc(B->size2);
	gsl_matrix_set_identity(B);
	for(int i=0; i<B->size2; i++){
		gsl_vector_view column = gsl_matrix_column(B,i);
		//put in a copy of these then solve
		gsl_blas_dcopy(&column.vector,x);
		GS_solve(Q,R,x,&column.vector);
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
	printf("\n");
	print_matrix("R=",R);
	printf("\n");
	print_matrix("Q=",Q);
	printf("\n");
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
	printf("We can show Qt*Q=1\n");
	print_matrix("Qt*Q=", QtimesQt);
	printf("\n");
	printf("I believe this is close to 0 and only isnt due to rounding errors\n");
	printf("\n");
	
	gsl_matrix * QR=gsl_matrix_alloc(n,m);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, Q,R,0.0, QR );
	printf("Checking QR=A\n");
	print_matrix("QR=",QR);
	printf("\n");

	printf("Now time for part 2 GS_solve\n");

	//put in the solve funct
	int nm=8;

	gsl_matrix* Ab=gsl_matrix_alloc(nm,nm);
	gsl_matrix* Ac=gsl_matrix_alloc(nm,nm);
	gsl_vector* b=gsl_vector_calloc(nm);
	gsl_vector* x=gsl_vector_calloc(nm);

	for(int i=0;i<Ab->size1;i++){
		for(int j=0;j<Ab->size2;j++){
			gsl_matrix_set(Ab,i,j,rand()%10);
		}
	}
	for(int i=0;i<b->size;i++){
		gsl_vector_set(b,i,rand()%10);
	}
	print_vector("b =",b);
	printf("\n");
	// factorise Ab into Qb and Rb
	
	gsl_matrix* Qb=gsl_matrix_alloc(nm,nm);
	gsl_matrix* Rb=gsl_matrix_alloc(nm,nm);
	gsl_matrix_memcpy(Qb, Ab);	

	GS_decomp(Qb,Rb);

	print_matrix("Qb=",Qb);
	printf("\n");
	print_matrix("Rb=",Rb);
	printf("\n");
	
	GS_solve(Qb,Rb,b,x);
	
	
	printf("showing Ax=b\n");
	//to print out b we need a print vector func
	gsl_vector * Ax=gsl_vector_calloc(nm);
	gsl_blas_dgemv(CblasNoTrans, 1, Ab,x,0.0, Ax );//this is the same as previous but vector 
	print_vector("A*x=",Ax);
	printf("\n");
	printf("to check against b\n");
	

	print_vector("b=",b);
	printf("As can be seen there has been a large amount of rounding and so Ax and b are not fully equal but they are close enough to show the function works as intended");
	printf("\n");
	printf("\n");
	printf("\n");
	printf("\n");

	printf("Part B\n");

	gsl_matrix* T=gsl_matrix_alloc(n,n);
	gsl_matrix* S=gsl_matrix_alloc(n,n);
	gsl_matrix* V=gsl_matrix_alloc(n,n);
	
	for(int i=0;i<T->size1;i++){
		for(int j=0;j<T->size2;j++){
			gsl_matrix_set(T,i,j,rand()%10);
		}
	}
	
	gsl_matrix_memcpy(S, T);
	GS_decomp(S,V);

	print_matrix("A=",T);
	print_matrix("Q=",S);
	print_matrix("R=",V);

	gsl_matrix *invT = gsl_matrix_alloc(n,n);

	GS_inverse(S,V,invT);
	print_matrix("Inverse of A =",invT);

	printf("To check we calculate AB = I = BA\n");
	
	gsl_matrix *AB = gsl_matrix_alloc(n,n);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, T,invT,0.0, AB);
	print_matrix("AB =",AB);

	gsl_matrix *BA = gsl_matrix_alloc(n,n);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, invT,T,0.0, BA);
	print_matrix("BA =",BA);





return 0;
}
