#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include<assert.h>
#include"leastsf.h"

double f(int i, double x){
	switch(i){
		case 0: return 1; break;
		case 1: return x ; break;
		case 2: return x*x; break;
		default: printf("f is wrong %i",i); return NAN; break;
	}
}



void LSF(gsl_vector * t, gsl_vector * y, gsl_vector * dy, int m, double (*f)(int i, double x), gsl_vector * c, gsl_matrix * S){
	int n = t->size;

	gsl_matrix* A = gsl_matrix_alloc(n,m);
	gsl_vector* b = gsl_vector_alloc(n);
	gsl_matrix* R = gsl_matrix_alloc(m,m);

	for(int i=0;i<n;i++){
		gsl_vector_set(b,i,gsl_vector_get(y,i)/gsl_vector_get(dy,i));
				
		for(int k=0;k<m;k++){
			gsl_matrix_set(A,i,k,f(k,gsl_vector_get(t,i))/gsl_vector_get(dy,i));
		}
	}
	GS_decomp(A,R);
//	GS_solve(A,R,b,c);//missed this step but you need to back subs and put into c like the booklet says	
	
	gsl_blas_dgemv(CblasTrans,1,A,b,0,c);
	for(int i=c->size-1;i>=0;i--){
		double s=gsl_vector_get(c,i);
		for(int k=i+1;k<c->size;k++){
			s-=gsl_matrix_get(R,i,k)*gsl_vector_get(c,k);
		}
		gsl_vector_set(c,i,s/gsl_matrix_get(R,i,i));
	}
	gsl_matrix* inverse = gsl_matrix_alloc(m,m);
	gsl_matrix* id = gsl_matrix_alloc(m,m);
	gsl_matrix_set_identity(id);//use from inverse previous problem
	GS_inverse(id,R,inverse);
	print_matrix("I =",id);
	print_matrix("R =",R);
	print_matrix("Inverse =",inverse);

	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, inverse,R,0.0, S);
	print_matrix("inverse*R=",S);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1.0, R,inverse,0.0, S);
	print_matrix("R*inverse=",S);

	//now wants to multiply inverse R and inverse transpose R
	gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1.0, inverse,inverse,0.0, S);

	gsl_matrix_free(A);
	gsl_matrix_free(R);
	gsl_matrix_free(inverse);
	gsl_vector_free(b);
	gsl_matrix_free(id);
}








//move the functions into  a .h as linear was horrible with that

int main(){

	//start by putting the data into vectors
	
	int n = 9;
	gsl_vector* t=gsl_vector_alloc(n);
	gsl_vector* y=gsl_vector_calloc(n);
	gsl_vector* dy=gsl_vector_calloc(n);
	
	//put the data in one at a time starting with t
	gsl_vector_set(t,0,1);
	gsl_vector_set(t,1,2);
	gsl_vector_set(t,2,3);
	gsl_vector_set(t,3,4);
	gsl_vector_set(t,4,6);
	gsl_vector_set(t,5,9);
	gsl_vector_set(t,6,10);
	gsl_vector_set(t,7,13);
	gsl_vector_set(t,8,15);

	//put in the y vector
	gsl_vector_set(y,0,117);
	gsl_vector_set(y,1,100);
	gsl_vector_set(y,2,88);
	gsl_vector_set(y,3,72);
	gsl_vector_set(y,4,53);
	gsl_vector_set(y,5,29.5);
	gsl_vector_set(y,6,25.2);
	gsl_vector_set(y,7,15.2);
	gsl_vector_set(y,8,11.1);

	//set the dy
	for(int i=0;i<9;i++){
		gsl_vector_set(dy,i,(gsl_vector_get(y,i)/20));
	}
	
	printf("Part A\n");
	printf("The vectors are:\n");
	print_vector("t=",t);
	print_vector("y=",y);
	print_vector("dy=",dy);
	//now its time to make the least square fit function


	int m=2;
	gsl_vector * c =gsl_vector_alloc(m);
	gsl_matrix *S = gsl_matrix_alloc(m,m);

	//switch to log 
	//uncertanty is dlny=dy/y
	gsl_vector_div(dy,y);//stores new dy n dy
	//gsl_vector_mul(y,log 
	for(int k=0; k<y->size; k++){
		gsl_vector_set(y,k,log((gsl_vector_get(y,k))));
	}
	print_vector("log y =",y);
	print_vector("dy =",dy);
	//now put into the least square func
	LSF(t,y,dy,m,f,c,S);
	
	print_matrix("S=",S);
	print_vector("c=",c);
	
	//the outfiles is very full so put into a new txt file to plot to check numbers are right
	FILE* fit = fopen("LSF.out.txt","w");
	FILE* plot = fopen("plot.out.txt","w");
	int p = 9;
	for(int i=0;i<p;i++){
		fprintf(plot,"%10g %10g %10g\n",gsl_vector_get(t ,i),gsl_vector_get(y,i), gsl_vector_get(dy,i));
	}
	int o = 100;
	gsl_vector* mt=gsl_vector_alloc(o);
	gsl_vector* my=gsl_vector_alloc(o);

	for(int i=0;i<o-1; i++){
		gsl_vector_set(mt,i,16*i/(o-1));
		gsl_vector_set(my,i,gsl_vector_get(c,0)+gsl_vector_get(mt,i)*gsl_vector_get(c,1));

		fprintf(fit, "%10g %10g\n",gsl_vector_get(mt,i), gsl_vector_get(my,i));
	}





	printf("Halflife of ThX measured is %g +- %g\n The modern value of 224Ra is 3.6 days and so we are slightly over that value. However, There are other Ra isotopes with longer halflives which may be presnt in the sample they used which increases the time measured..",-log(2)/gsl_vector_get(c,1),log(2)/pow(gsl_vector_get(c,1),2)*sqrt(gsl_matrix_get(S,1,1)));

	printf("\n\n\n\n\n");
//If I come back and dont know what I did, I couldnt get this to work but it turns out I had coppied a section of my code wrng and this is a check 
//
	
//	gsl_matrix* A=gsl_matrix_alloc(m,m);
//	for(int i=0;i<A->size1;i++){
//		for(int j=0;j<A->size2;j++){
//			gsl_matrix_set(A,i,j,rand()%10);
//		}
//	}

//	print_matrix("A=",A);
	
//	gsl_matrix* inv = gsl_matrix_alloc(m,m);
//	gsl_matrix* ide = gsl_matrix_alloc(m,m);
//	gsl_matrix_set_identity(ide);
//	GS_inverse(ide,A,inv);
//
//	print_matrix("inv = ",inv);





return 0;
}
