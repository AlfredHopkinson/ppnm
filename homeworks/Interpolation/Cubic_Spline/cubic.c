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



typedef struct {int n; double *x, *y, *b, *c, *d;} cubicspline;

cubicspline * cubicspline_alloc(int n, double *x, double *y){ 
	cubicspline *s = (cubicspline*)malloc(sizeof(cubicspline));
	s->x = (double*) malloc(n*sizeof(double));
	s->y = (double*) malloc((n*sizeof(double));
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
	D[0]=2 for(int i=0;i<n-2;i++)D[i+1] = 2*h[i]/h[i+1]+2; D[n-1]=2;
	Q[0]=1 for(int i=0;i<n-2;i++) 
	
	
	
	s->c[0]=0;
	for (i=0;i<n-2;i++) {s->c[i+1] = (p[i+1]-p[i]-s->c[i]*h[i])/h[i+1];}
	s->c[n-2]/=2;
	for (i=n-3;i>=0;i--) {s->c[i] = (p[i+1]-p[i]- s->c[i+1]*h[i+1])/h[i];}
	for(i=0;i<n-1;i++) {s->b[i]=p[i]-s->c[i]*h[i];}
	return s;
}


