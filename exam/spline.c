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
#include<float.h>
#include"spline.h"

//include the binsearch from previous problem
int binsearch(int n, double* x, double z){ 
	assert(n>1 && z>=x[0] && z<=x[n-1]);
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
}





//exam says take insperation from akima seen as that is a subspline interpolation. So I am using that as a basic structure

typedef struct {int n; double *x,*y,*c,*d;} akima_spline;
akima_spline* akima_spline_alloc(int n, double *x, double *y){
	assert(n>2); double h[n-1],p[n-1];
	for(int i=0;i<n-1;i++){h[i]=x[i+1]-x[i]; assert(h[i]>0);}
	for(int i=0;i<n-1;i++){p[i]=(y[i+1]-y[i])/h[i];}
	akima_spline *s = (akima_spline*)malloc(sizeof(akima_spline));
	s->x = (double*)malloc(n*sizeof(double));
	s->y = (double*)malloc(n*sizeof(double));
	s->b = (double*)malloc(n*sizeof(double));
	s->c = (double*)malloc((n-1)*sizeof(double));
	s->d = (double*)malloc((n-1)*sizeof(double));
	s->n =n; for(int i=0;i<n;i++){s->x[i]=x[i]; s->y[i]=y[i];}
	s->b[0] = p[0]; s->b[1] = (p[0]+p[1])/2;
	for(int i=2;i<n-2;i++){
		double w1=fabs (p[i+1]-p[i]), w2=fabs(p[i-1]-p[i-2]);
		if (w1+w2==0) s->b[i]=(p[i+1]+p[i])/2;
		else s->b[i]=(w1*p[i-1]+w2*p[i])/(w1+w2);
	}
	for (int i=0;i<n-1;i++){
		s->c[i]=(3*p[i]-2*s->b[i]-s->b[i+1])/h[i];
		s->d[i]=(s->b[i+1]+s->b[i]-2*p[i])/h[i]/h[i];
	}
	return s;
}

double akima_spline_eval(akima_spline *s, double z){
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	int i = 0, j=s->n-1;
	while (j-i>1){int m=(i+j)/2; if (z>s->x[m]) i=m; else j=m;}
	double h=z-s->x[i];
	return s->y[i]+h*(s->b[i]+h*(s->c[i]+h*s->d[i]));
}

void akima_spline_free(akima_spline *s){
	free(s->x);
	free(s->y);
	free(s->b);
	free(s->c);
	free(s->d);
	free(s);
}
