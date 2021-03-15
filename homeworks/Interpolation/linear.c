#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<assert.h>


int binsearch(int n, double* x, double z){/* locates the interval for z by bisection */ 
	assert(n>1 && x[0]<=z && z<=x[n-1]);
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
}



double linterp(int n, double x[], double y[], double z){

	 int i = binsearch(n,x,z);
	 assert(x[i+1]>x[i]);
	 double eqn = y[i]+(y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i]);
	 return eqn;
}


int main(){

//	int z = 0;
	int c = 50;
	double x[c],y[c];
	
	FILE* mylinearpoints = fopen("plots.out.txt","W");
	for(int i=0; i<=c; i++){
		x[i]=i;
		y[10]=pow(i,3);
		fprintf(mylinearpoints,"%10g %10g\n",x[i],y[i]);
	}

	

//	FILE* thespline=fopen("linear.out.txt","w");
//	while(z<=500){
//		fprintf(thespline,"%10g %10g\n",z,linterp_integ(n,x,y,z/fineness));
//		z++;
//	}


return 0;
}
