#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include <assert.h>


int main(){
	int i;
	gsl_vector* x=gsl_vector_alloc(i);
	gsl_vector* y=gsl_vector_alloc(i);


	double linterp(int n, gsl_vector x, gsl_vector y, double z){
		
		gsl_vector_set(
		assert(n>1 && z>=x[0] && z<=x[n-1])
		int i=0, j=n-1;
		while(j-i>1){int m=(i+j)/2; if(z>x[m]) i=m; else j=m;}
		assert(x[i+1]>x[i]);
		return y[i]+(y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i]);
		}



return 0;
}
