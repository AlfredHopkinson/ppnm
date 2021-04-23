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

void GS_decomp(gsl_matrix *A, gsl_matrix *R){
	int m = A->size2;
	int n = A->size1;
	assert(A->size2 == R ->size1); 
	double Rij, Rmatrixii;
	for(int i=0;i<m;i++){
		gsl_vector_view icolumn = gsl_matrix_column(A, i);
		Rmatrixii= gsl_blas_dnrm2(&icolumn.vector); 
		gsl_matrix_set(R, i, i, Rmatrixii);
		gsl_vector_scale(&icolumn.vector,1.0/Rmatrixii);
				
		for(int j=i+1;j<m;j++){
			gsl_vector_view columnj = gsl_matrix_column(A,j);
			gsl_blas_ddot(&icolumn.vector, &columnj.vector, &Rij);	 
			gsl_matrix_set(R,i,j,Rij);
			gsl_blas_daxpy(-Rij, &icolumn.vector, &columnj.vector);
		}
	}
}




void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x) {
	gsl_blas_dgemv(CblasTrans,1,Q,b,0,x);
	for(int i=x->size-1; i>=0; i--) {
		double s=gsl_vector_get(x,i);
		for(int k=i+1; k<x->size; k++) s-=gsl_matrix_get(R,i,k)*gsl_vector_get(x,k);
			gsl_vector_set(x,i,s/gsl_matrix_get(R,i,i));
	}
}

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
	gsl_vector*x = gsl_vector_alloc(B->size2);
	gsl_matrix_set_identity(B);
	for(int i=0; i<B->size2; i++){
		gsl_vector_view column = gsl_matrix_column(B,i);
		gsl_blas_dcopy(&column.vector,x);
		GS_solve(Q,R,x,&column.vector);
	}
}





//void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
//	gsl_vector*x = gsl_vector_alloc(B->size2);
//	gsl_matrix_set_identity(B);
//	for(int i=0; i<B->size2; i++){
//		gsl_vector_view column = gsl_matrix_column(B,i);
//		gsl_blas_dcopy(&column.vector,x);
//		GS_solve(Q,R,x,&column.vector);
//	}
//}

// Moved al;l this to main as it wasnt working in here for some reason

//here we were given a python least square fit so convert it to c or use basis to make own

//void LSF(gsl_vector * t, gsl_vector * y, gsl_vector * dy, int m, double (*f)(int i, double x), gsl_vector * c, gsl_matrix * S){
//	int n = t->size;

//	gsl_matrix* A = gsl_matrix_alloc(n,m);
//	gsl_vector* b = gsl_vector_alloc(n);
//	gsl_matrix* R = gsl_matrix_alloc(m,m);

//	for(int i=0;i<n;i++){
//		gsl_vector_set(b,i,gsl_vector_get(y,i)/gsl_vector_get(dy,i));
		
//		for(int k=0;k<m;k++){
//			gsl_matrix_set(A,i,k,f(k,gsl_vector_get(t,i))/gsl_vector_get(dy,i));
//		}
//	}
//	GS_decomp(A,R);
//	GS_solve(A,R,b,c);//missed this step but you need to back subs and put into c like the booklet says
	
//	gsl_blas_dgemv(CblasTrans,1,A,b,0,c);
//	for(int i=c->size-1;i>=0;i--){
//		double s=gsl_vector_get(c,i);
//		for(int k=i+1;k<c->size;k++){
//			s-=gsl_matrix_get(R,i,k)*gsl_vector_get(c,k);
//		}
//		gsl_vector_set(c,i,s/gsl_matrix_get(R,i,i));
		
//	}

//	gsl_matrix* inverse = gsl_matrix_alloc(m,m);

//	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, inverse,R,0.0, S);
//	print_matrix("inverse*R=",S);
//	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, R,inverse,0.0, S);
//	print_matrix("R*inverse=",S);


//	gsl_matrix* id = gsl_matrix_alloc(m,m);
//	gsl_matrix_set_identity(id);//use from inverse previous problem
	
	//GS_inverse(id,R,inverse);
	//now wants to multiply inverse R and inverse transpose R
//	gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, inverse,inverse,0.0, S);

//	gsl_matrix_free(A);
//	gsl_matrix_free(R);
//	gsl_matrix_free(inverse);
//	gsl_vector_free(b);
//	gsl_matrix_free(id);	

//}




