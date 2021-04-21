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

double integrate(double (*f)(double), double a, double b, double d, double e, double f2, double f3, int nrec)
{

	double f1 = f(a+(b-a)/6);
	double f4 = f(a+5*(b-a)/6);


	double Q = Q=(2*f1+f2+f3+2*f4)/6*(b−a);
	double q=(f1+f4+f2+f3)/4*(b−a);
       	double tol = d+e*fabs(Q);
	double err= fabs(Q-q); //not sure about the rules here but going to use the ones from the example code
	if (err < d+e*fabs(Q)) return Q;
	else return integrate(f,a,(a+b)/2,d/√2,e)+ integrate(f,(a+b)/2,b,d/√2,e);
}

double recadapt(double (*f)(double), double a, double b, double d, double e){
	double f2 = f(a+2*(b-a)/6);
	double f3 = f(a+4*(b-a)/6);
	int nrec = 0;
	return integrate(f,a,b,d,e,f2,f3,nrec);
}

