#include <math.h>
#include <complex.h>
#include <stdio.h>

int main(){
	double gamm = tgamma(5);
	double bess = j1(0.5);

	complex c = csqrt(-2);
	complex d = cexp(I*M_PI);
	complex e = cexp(I);
	complex f = cpow(I,M_E);
	complex g = cpow(I,I);
	
	float x_float = 1.f/9;
	double x_double = 1./9;
	long double x_long_double = 1.L/9;
	
	printf("gamma(5) = %g\n", gamm);
	printf("bessel(0.5) = %g\n", bess);

	printf("sqrt(-2) = %g + I%g\n", creal(c), cimag(c));
	printf("exp(iPI) = %g + I%g\n", creal(d), cimag(d));
	printf("exp(i) =  %g + I%g\n", creal(e), cimag(e));
	printf("i^e =  %g + I%g\n", creal(f), cimag(f));
	printf("i^i =  %g + I%g\n", creal(g), cimag(g));

	printf("float output for 1/9 = %.25g", x_float);
	printf("double output for 1/9 = %.25lg", x_double);
	printf("Long double output for 1/9 = %.25Lg", x_long_double);

return 0;
}

