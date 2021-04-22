#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include<assert.h>
#include <gsl/gsl_integration.h>
//#include int.h I cant get this to work it rerads a void 
//just going back to put it in main
//thinking about it it doesnt like those functions in a function so that was probably it


double integrate(double (*f)(double), double a, double b, double d, double e, double f2, double f3, int nrec)
{

	double f1 = f(a+(b-a)/6);
	double f4 = f(a+5*(b-a)/6);

	double Q = (b-a)/6*(2*f1+f2+f3+2*f4);

	double q = (b-a)/4*(f1+f2+f3+f4);
       	double tol = d+e*fabs(Q);
	double err= fabs(Q-q); //not sure about the rules here but going to use the ones from the example code
	if (err < d+e*fabs(Q)) return Q;
		else return integrate(f,a,(a+b)/2,f1, f2,d/sqrt(2),e,nrec+1)+ integrate(f,(a+b)/2,b,f3,f4,d/sqrt(2),e,nrec+1);
}



double recadapt(double (*f)(double), double a, double b, double d, double e){
	double f2 = f(a+2*(b-a)/6);
	double f3 = f(a+4*(b-a)/6);
	int nrec = 0;
	return integrate(f,a,b,d,e,f2,f3,nrec);
}



//here choosing the second integral given as an example 
double f(double x){return sqrt(x);}

int main(){
	printf("Part A\n");


	double a=0.0;
	double b = 1.0;
	double d = 0.01; //not sure as to my accuracy but went for somthing reasonable
	double e = 0.01;
	int nrev = 0;

	double function(double x){nrev++; return f(x);}
	double integ = recadapt(function,a,b,d,e);

	printf("The Test Integrals\n");
	printf("Integral of f(x)=sqrt(x) is %g and it used %i number of evaluations. It correct value is 2/3 or 0.6666 so it is very close.\n",integ,nrev);


return 0;
}
