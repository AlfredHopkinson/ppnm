#include<stdio.h>
#include<limits.h>
#include<stdlib.h>
#include<assert.h>
#include<float.h>
#include<math.h>

int equal(double a, double b, double tau, double epsilon){
	if(fabs(a-b)< tau || (fabs(a-b))/(fabs(a)+fabs(b)) < tau/2) return 1;
return 0;
}
