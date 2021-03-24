#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include<assert.h>
#include"leastsf.h"

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

return 0;
}
