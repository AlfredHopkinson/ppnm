#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include<assert.h>
#include"func.h"

//example that is in the exercise is a sho
void example(double t, gsl_vector* y, gsl_vector* dydt){
	gsl_vector_set(dydt,0,gsl_vector_get(y,1));
	gsl_vector_set(dydt,1,-gsl_vector_get(y,0));
}



int main(){
	//firslty make the vectors 
	gsl_vector* yt = gsl_vector_alloc(3);
	gsl_vector* yh = gsl_vector_alloc(3);
	gsl_vector* err = gsl_vector_alloc(3);
	//set the initial conditions
	gsl_vector_set(yt,0,1);
	gsl_vector_set(yt,1,0);
	
	//now use the driver for d2u=-u
	double h =0.5;

	driver(example,0,yt,4*M_PI,yh,err,h,0.1,0.1);

return 0;
}
