#include<stdio.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#include<assert.h>
#include<gsl/gsl_interp.h>
#include<gsl/gsl_spline.h>
#include<stdlib.h>


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

	
	double integral = 0;
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid;
		else j=mid;
	}
	for(int k=1;k<=i;k++){ 
		integral+=y[k-1]*(x[k]-x[k-1])+0.5*(y[k]-y[k-1])*(x[k]-x[k-1]);
		}	

	double yin = linterp(n,x,y,z);
	integral +=y[i]*(z-x[i])+0.5*(yin-y[i])*(z-x[i]);

	return integral;
}


double linterp_gsl(int n, double* x, double* y, double z){
	gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_linear,n);
	gsl_interp_init(interpolation,x,y,n);
	gsl_interp_accel*accelerator = gsl_interp_accel_alloc();

return gsl_interp_eval(interpolation,x,y,z,accelerator);
}

double linterp_integ_gsl(int n, double* x, double* y, double z){
int i=0, j=n-1;
while(j-i>1){
	int mid=(i+j)/2;
	if(z>x[mid]) i=mid;
	else j=mid;
	}
gsl_interp *interpolation = gsl_interp_alloc(gsl_interp_linear,n);
gsl_interp_init(interpolation,x,y,n);
gsl_interp_accel*accelerator = gsl_interp_accel_alloc();
double integral=0;
integral+= gsl_interp_eval_integ(interpolation,x,y,0,z,accelerator);

gsl_interp_free(interpolation);
return  integral;
}





//this didnt work so try simplifying it and using a while loop
//	for(p=1; p<=i; p++){
//		if(x[i+1]>x[i]){
//			double square = (x[i+1]-x[i])*y[i];
//			double triangle = 0.5*(x[i+1]-x[i])*(y[i+1]-y[i]);
//			double partintegral = square + triangle;
//			integral += partintegral;
//		}
//		else{
//			double square = (x[i+1]-x[i])*y[i+1];
//			double triangle = 0.5*(x[i+1]-x[i])*(y[i]-y[i+1]);
//			double partintegral = square + triangle;
//			integral += partintegral;
//		}
//		}
//	return integral;





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
	fclose(mylinearpoints);

//now i need to interpolate which we will do putting in smaller values through the function
	double thin = 5;
	while(z<=thin*x[n-1]){
		printf("%10g %10g\n",z/thin,linterp(n,x,y,z/thin));
		z++;
	}

	//to integrate we want to put it into a new file 
	FILE* theinteg=fopen("linearint.out.txt","w");
	for(z=0; z<=thin*(n-1);z++){
		fprintf(theinteg,"%10g %10g\n",z/thin,linterp_integ(n,x,y,z/thin));
	}

	//now I have to compare what I have done above with GSL own
	
	


	FILE* gsldata = fopen("gslinfo.txt","w");	
	for(z=0; z<=thin*(n-1);z++){
		fprintf(gsldata,"%10g %10g %10g \n",z/thin, linterp_gsl(n,x,y,z/thin),linterp_integ_gsl(n,x,y,z/thin));
	}



return 0;
}
