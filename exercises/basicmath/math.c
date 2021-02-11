#include <stdio.h>
#include <math.h>
#include <complex.h>

int main(){
	double gamm = tgamma(5);
	double bess = j1(0.5);

	complex root2 = csqrt(-2);
	complex expipi = cexp(I*M_PI);
	complex expi = cexp(I);
	complex ipowere = cpow(I,M_E);
	complex ipoweri = cpow(I,I);

	float x_float = 1.f/9;
	double x_double = 1./9;
	long double x_long_double = 1.L/9;

	printf("gamma(5) = %g\n", gamm);
	printf("bessel(0.5) = %g\n", bess);

	printf("sqrt(-2) = %g+%g\n", creal(root2), cimag(root2));
	printf("exp((iPI) = %g+%g\n", creal(expipi), cimag(expipi));
	printf("exp(i) = %g+%g\n", creal(expi), cimag(expi));
	printf("i^e = %g+%g\n", creal(ipowere), cimag(ipowere));
	printf("i^i = %g+%g\n", creal(ipoweri), cimag(ipoweri));

	printf("1/9 as a float = %.25g", x_float);
	printf("1/9 as a double = %.25lg", x_double);
	printf("1/9 as a long double = %.25Lg", x_long_double);

}
