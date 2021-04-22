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

//take from prev
//get the rest from the webpage then the derivitive and int can be put in after
//then add the gsl sub routines and plot them over this one



typedef struct {int n; double *x, *y, *b, *c, *d;} cubic_spline;

cubic_spline * cubic_spline_alloc(int n, double *x, double *y){ 
	cubic_spline *s = (cubic_spline*)malloc(sizeof(cubic_spline));
	s->x = (double*) malloc(n*sizeof(double));
	s->y = (double*) malloc(n*sizeof(double));
	s->b = (double*) malloc(n*sizeof(double));
	s->c = (double*) malloc((n-1)*sizeof(double));
	s->d = (double*) malloc((n-1)*sizeof(double));
	s->n = n;
	int i;
	for(i=0;i<n;i++){s->x[i]=x[i];s->y[i]=y[i];}
	double p[n-1], h[n-1];
	for(i=0;i<n-1;i++){h[i]=x[i+1]-x[i]; assert(h[i]>0);}
	for(i=0;i<n-1;i++)p[i]=(y[i+1]-y[i])/h[i];
	double D[n], Q[n-1], B[n];
	D[0]=2;
	Q[0]=1;	
	for(int i=0;i<n-2;i++){D[i+1] = 2*h[i]/h[i+1]+2; D[n-1]=2;}
	for(int i=0;i<n-2;i++){Q[i+1] = h[i]/h[i+1];} 
	for(int i=0;i<n-2;i++){B[i+1] = 3*(p[i]+p[i+1]*h[i]/h[i+1]);}
	B[0] = 3*p[0];
	B[n-1]=3*p[n-2];
	for(int i=1;i<n;i++){D[i]-=Q[i-1]/D[i-1];B[i]-=B[i-1]/D[i-1];}
	s->b[n-1]=B[n-1]/D[n-1];
	for(int i=n-2;i>=0;i--){ s->b[i]=(B[i]-Q[i]*s->b[i+1])/D[i];}
	for(int i=0;i<n-1;i++){
		s->c[i]=(-2*s->b[i]-s->b[i+1]+3*p[i])/h[i];
		s->d[i]=(s->b[i]+s->b[i+1]-2*p[i])/h[i]/h[i];
	}
return s;
}

	
	


double cubic_spline_eval(cubic_spline *s, double z){
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	int i=0, j=s->n-1;
	while(j-i>1){int m=(i+j)/2;
		if (z>s->x[m]) i=m;
		else j=m;}
	double h=z-s->x[i];
	return s->y[i]+h*(s->b[i]+h*(s->c[i]+h*s->d[i]));
}

//now basically use the int and deriv basis from the quad here
double cspline_derivitive(cubic_spline *s, double z){
	assert (z>=s->x[0] && z<=s->x[s->n-1]); 
	int i = binsearch(s->n,s->x,z);
	double dx = z - s->x[i];
	return s->b[i] + 2*s->c[i]*dx;
}

double cspline_integ(cubic_spline *s, double z){
	assert (z>=s->x[0] && z<=s->x[s->n-1]);
	int i = binsearch(s->n,s->x,z);
	double integ;
	double dx;
	for (int j = 0; j < i; j++){
		dx = s->x[j+1] - s->x[j];
		integ += dx*(s->y[j] + dx*(s->b[j]/2 + dx*s->c[j]/3));
	}
	dx=z-s->x[i];
	integ += dx*(s->y[i]+dx*(s->b[i]/2+dx*s->c[i]/3));
	return integ;
}



void cubic_spline_free (cubic_spline *s){
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s->d);
	free(s);
}






int main(){
	//need the allocation of the points 
	int n=10;
	double x[n], y[n];
	FILE* mypoints = fopen("plots.out.txt","w");
	
	for(int i=0;i<n;i++){
		x[i]=i;
		y[i]= (100.0*rand()/RAND_MAX);
		fprintf(mypoints,"%10g %10g\n",x[i],y[i]);
	}
	
	double thin = 10;
	int z = 0; //taken and modified from linear
	cubic_spline * s = cubic_spline_alloc(n,x,y);
		
	FILE* cubout=fopen("cub.out.txt","w");
	for(z=0; z<=thin*(n-1);z++){
		fprintf(cubout,"%10g %10g %10g %10g\n",z/thin, cubic_spline_eval(s,z/thin), cspline_derivitive(s,z/thin), cspline_integ(s,z/thin));
	}
return 0;		
}







