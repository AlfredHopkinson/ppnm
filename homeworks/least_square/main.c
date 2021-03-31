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
	FILE* plot = fopen("LSF.out.txt","w");
	int o = 100;
	gsl_vector* mt=gsl_vector_alloc(o);
	gsl_vector* my=gsl_vector_alloc(o);

	for(int i=0;i<o-1; i++){
		gsl_vector_set(mt,i,16*i/(o-1));
		gsl_vector_set(my,i,gsl_vector_get(c,0)+gsl_vector_get(mt,i)*gsl_vector_get(c,1));

		fprintf(plot, "%10g %10g\n",gsl_vector_get(mt,i), gsl_vector_get(my,i));
	}





	printf("halflife of ThX is %g +- %g\n",-log(2)/gsl_vector_get(c,1),log(2)/pow(gsl_vector_get(c,1),2)*sqrt(gsl_matrix_get(S,1,1)));

return 0;
}
