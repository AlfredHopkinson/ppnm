#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include<assert.h>
#include <gsl/gsl_integration.h>
#define RND ((double)rand()/RAND_MAX)


// I firstly tried to do this with the example on the question but I will actually now go back and restart using the chapter for help as its broken down into parts a lot more
//void plainmc(int dim,double f(int dim,double* x),double* a,double* b,int N){
//	double V=1; for(int i=0;i<dim;i++)V*=b[i]-a[i];
//	double sum=0,sum2=0,x[dim];
//	for(int i=0;i<N;i++){
//		for(int i=0;i<dim;i++)x[i]=a[i]+RANDOM*(b[i]-a[i]);
//			double fx=f(dim,x); sum+=fx; sum2+=fx*fx;
//	}
//	double mean=sum/N, sigma=sqrt(sum2/N-mean*mean);
//	result= mean*V+sigma*V/sqrt(N);
//return result;
//}

//Make a random first template in chapter
void rand(int d, double *a, double *b, double *x){
	for(int i=0;i<d;i++){
		x[i]=a[i]+RND*(b[i]-a[i]);
	}
}

//using the chapter template this time
void plainmc(int dim,double f(int dim,double* x),double* a,double* b,int N, double*result, double*error){
	double V=1;
	for (int i=0;i<dim;i++){
		V*=b[i]-a[i];
	}
	double sum=0, sum2=0,fx,x[dim];
	for(int i=0;i<N;i++){
		rand(dim,a,b,x);
		fx=f(dim,x);
		sum+=fx;
		sum2+=fx*fx;
	}
	double avr =sum/N, var = sum2/N-avr*avr;
	*result = avr*V;
	*error = sqrt(var/N)*V;
}


int main(){
	//need to do some test integrals 
	//try to int a simple square
	printf("Example integral of a square with sides of a length 2,4\n");
	double res, err;
	
	double a[]={0, 0}, b[]={2,4};
	double square;
	double f(int dimens, double *x){
		
		square = x[0];
		return square;
	}
	int dimens = 2, N = 1e3;
	double squarearea = 2*4;
	plainmc(dimens,f,a,b,N,&res,&err);
	printf("The real area of the square is %g the calculated area is %g, the associated error is %g\n",squarearea,res,err);

	printf("\n");
	printf("\n");

	//do a sphere of radius 1
	printf("Example Integral of a sphere of radius 1\n");
	double result, error;
	double a1[]={-1, -1, -1};
	double b1[]={1,1,1};//set the interation values
	double g(int dimens, double * x) {
		double value = 1.0, sum = 0.0, norm; 
		for (int i = 0;i<dimens;i++) sum+=x[i]*x[i]; 
		norm = sqrt(sum); 
		if (norm>1) value*=0.0; //set g
	return value;
	}
	
	int dimens1 = 3, N1 = 1e6;
	plainmc(dimens1,g,a1,b1,N1,&result,&error);
	printf("Real volume = %g, calculated volume = %g, associated error = %g\n",4*M_PI*pow(1,3)/3,result,error);
	
	printf("\n");
	printf("\n");


	printf("Integral from the homework\n");
	double result2, error2;
	double a2[]={0,0,0}, b2[]={M_PI,M_PI,M_PI};
	double h(int dim, double *x){
		double func = 1-cos(x[0])*cos(x[1])*cos(x[2]);
		return pow(func,-1)*(1/pow(M_PI,3));
	}
	int dimens2 = 3, N2 = 1e6;
	plainmc(dimens2,h,a2,b2,N2,&result2,&error2);

	double manual = 1.3932039296856768591842462603255;

	printf("Integral value = %g, Calculated value = %g, associated error = %g\n",manual, result2, error2);

	printf("\n");
	printf("\n");












return 0;
}
