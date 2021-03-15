#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<assert.h>


int binsearch(int n, double* x, double z){/* locates the interval for z by bisection */ 
	assert(n>1 && z>=x[0] && z<=x[n-1]);
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
}


//this is the linear interpolation with the intervals from binsearch
double linterp(int n, double* x, double* y, double z){

	 int i = binsearch(n,x,z);
	 assert(x[i+1]>x[i]);
	 double eqn = y[i]+(y[i+1]-y[i])/(x[i+1]-x[i])*(z-x[i]);
	 return eqn;
}

//analytiucal integration means it will be the area of the square and then the triangle
double linterp_integ(int n, double* x, double* y, double z){

	int i = binsearch(n,x,z);
	double area=0;
	int c=0;

	while(c<=i){
		area+=y[c-1]*(x[c]-x[c-1])+0.5*(y[c]-y[c-1])*(x[c]-x[c-1]);
		c++;
	}

	double yin = linterp(n,x,y,z);
	area+=y[i]*(z-x[i])+0.5*(yin-y[i])*(z-x[i]);

	return area;
}



int main(){
	
	int n=10;
	int z = 0;
//	here plot points first, started trying to manually put them in but decided to use the randdom function from previous excercises

	double x[n],y[n];
//make a plot of these points here when i can
	FILE* mylinearpoints = fopen("plots.out.txt","w");
	for(int i=0; i<n; i++){
		x[i]=i;
//		y[i]= (100.0*rand()/RAND_MAX);  //adapted from multi proc
		y[i]=pow(i,2);
		fprintf(mylinearpoints,"%10g %10g\n",x[i],y[i]);
	}

	double thin = 10;

	while(z<=thin*x[n-1]){
		printf("%10g %10g\n",z/thin,linterp(n,x,y,z/thin));
		z++;
	}

	FILE* theinteg=fopen("linearint.out.txt","w");
	while(z<=thin*x[n-1]){
		fprintf(theinteg,"%10g %10g\n",z/thin,linterp_integ(n,x,y,z/thin));
		z++;
	}



return 0;
}
