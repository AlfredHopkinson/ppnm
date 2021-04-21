#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include<assert.h>
#include <gsl/gsl_integration.h>
#include"int.h"

//here choosing the second integral given as an example 
double f(double x){return 4*sqrt(1-pow(x,2));}

int main(){
	printf("Part A\n");


	double a=0.0;
	double b = 1.0;
	double d = 0.001; //not sure as to my accuracy but went for somthing reasonable
	double e = 0.001;
	int nrev = 0;

	double function(double x){nrev++; return f(x);}
	double integ = recadapt(function,a,b,d,e);

	printf("The Test Integrals\n");
	printf("Integral of f(x)=4*sqrt(1-x^2) is %g and it used %i number of evaluations\n",integ,nrev);


return 0;
}
