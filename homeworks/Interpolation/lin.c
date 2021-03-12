#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>

int binsearch(int n, gsl_vector *x, double z);

double linterp(gsl_vector x, gsl_vector y, double z){
	int n = (x->size);
	int binsearch(n,x,z);

	double xi=gsl_vector_get(x,i);
	double xiplusone=gsl_vector_get(x,i+1);
	double yi=gsl_vectttor_get(y,i);
	double yiplusone=gsl_vector_get(y,i+1);

	double interp = yi + ((yiplusone-yi)/(xiplusone-xi))*(z-xi);

	return interp;
}

