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
//#include"spline.h"


typedef struct {int n; double *x,*y,*b,*c,*d;} akima_spline;
akima_spline* akima_spline_alloc(int n, double *x, double *y){
	assert(n>2); double h[n-1],p[n-1];
	for(int i=0;i<n-1;i++){h[i]=x[i+1]-x[i]; }
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
//	assert(z>=s->x[0] && z<=s->x[s->n-1]);
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
/*
//now I will try to create a cubic subspline that takes inspiration from this
typedef struct {int n; double *x,*y,*dy,*b,*c,*d;} cubic_spline;
cubic_spline* cubic_spline_alloc(int n, double *x, double *y, double *dy){

//	assert(n>2); 
	double dx[n-1],p[n-1],dp[n-1];
	for(int i=0;i<n-1;i++){
		dx[i] = x[i+1]-x[i];
	//	assert(dx[i]>0);
		p[i] = (y[i+1]-y[i])/dx[i];
		dp[i] = (dy[i+1]-dy[i])/dx[i];//gradient
	}
	
	cubic_spline *cs = (cubic_spline*)malloc(sizeof(cubic_spline));
	cs->x = (double*)malloc((n)*sizeof(double));
	cs->y = (double*)malloc((n)*sizeof(double));
	cs->b = (double*)malloc((n-1)*sizeof(double));
	cs->c = (double*)malloc((n-1)*sizeof(double));
	cs->d = (double*)malloc((n-1)*sizeof(double));
	


//	cs->n =n; for(int i=0;i<n-1;i++){cs->x[i]=x[i]; cs->y[i]=y[i];}

//	cs->b[0] = dp[0]; cs->b[1] = (dp[0]+dp[1])/2;
//	for(int i=2;i<n-2;i++){
//		double w1=fabs (dp[i+1]-dp[i]), w2=fabs(dp[i-1]-dp[i-2]);
//		if (w1+w2==0) cs->b[i]=(dp[i+1]+dp[i])/2;
//		else cs->b[i]=(w1*dp[i-1]+w2*dp[i])/(w1+w2);
//	}

	//b will have the gradient but the last point isnt needed
	for(int i=0; i<n-1; i++){cs->b[i]=dy[i];}

	//now the c and d coeficients **** check I got these right*****
	
	for(int i=0; i<n-1; i++){
		cs->c[i] = (3*(p[i]-cs->b[i])/dx[i])-dp[i];
		cs->d[i] = (dp[i]-2*cs->c[i])/(3*dx[i]);
	}
	return cs;
}
*/
typedef struct {int n; double *x,*y,*dy,*b,*c,*d;} cubic_spline;
cubic_spline* cubic_spline_alloc(int n, double *x, double *y, double *dy){
	double h[n-1],p[n-1],dp[n-1];
	for(int i=0;i<n-1;i++){h[i]=x[i+1]-x[i]; }
	for(int i=0;i<n-1;i++){p[i]=(y[i+1]-y[i])/h[i];}
	for(int i=0;i<n-1;i++){dp[i] = (dy[i+1]-dy[i])/h[i];}


	cubic_spline *cs = (cubic_spline*)malloc(sizeof(cubic_spline));
	
	cs->x = (double*)malloc(n*sizeof(double));
	cs->y = (double*)malloc(n*sizeof(double));
	cs->dy = (double*)malloc(n*sizeof(double));
	cs->b = (double*)malloc((n-1)*sizeof(double));
	cs->c = (double*)malloc((n-1)*sizeof(double));
	cs->d = (double*)malloc((n-1)*sizeof(double));
	

	cs->n =n; for(int i=0;i<n;i++){cs->x[i]=x[i]; cs->y[i]=y[i]; cs->dy[i] = dy[i];}
	cs->b[0] = dp[0]; cs->b[1] = (dp[0]+dp[1])/2;
//	for(int i=2;i<n-2;i++){
//		double w1=fabs (dp[i+1]-dp[i]), w2=fabs(dp[i-1]-dp[i-2]);
//		if (w1+w2==0) cs->b[i]=(dp[i+1]+dp[i])/2;
//		else cs->b[i]=(w1*dp[i-1]+w2*dp[i])/(w1+w2);
//	}
	printf("part 2 sub");
	for (int i=0;i<n-1;i++){
		cs->c[i] = (3*(p[i]-cs->b[i])/h[i])-dp[i];
		cs->d[i] = (dp[i]-2*cs->c[i])/(3*h[i]);
	}
	return cs;
}





double cubic_spline_eval(cubic_spline *cs, double z){
	assert(z>=cs->x[0] && z<=cs->x[cs->n-1]);
	int i = 0, j=cs->n-1;
	while (j-i>1){int m=(i+j)/2; if (z>cs->x[m]) i=m; else j=m;}
	double h=z-cs->x[i];
	return cs->y[i]+h*(cs->b[i]+h*(cs->c[i]+h*cs->d[i]));
}


void cubic_spline_free(cubic_spline *cs){
	free(cs->x);
	free(cs->y);
	free(cs->b);
	free(cs->c);
	free(cs->d);
	free(cs);
}









int main(){
//made some example data for the data set and I wanted it to be like the example in the chapter
//	double x[] = {0,1,2,3,4,5,6,7,8};
	
	double y[] = {-1,-1,-1,0,1,1,1};
	double dy[] = {0,0,0,1,0,0,0}; //a simple derivitive to start
	int n = (7);
	double x[n];
	printf("part 1");
	FILE* simplepoints = fopen("simple_points.txt","w");
	for (int i=0; i<n; i++){
		x[i] = i;
		fprintf(simplepoints, "%10g %10g %10g\n", x[i], y[i], dy[i]);
	}
	


	akima_spline *s = akima_spline_alloc(n,x,y);

	cubic_spline *cs = cubic_spline_alloc(n,x,y,dy);
	
	FILE* simple_out = fopen("simple_out.txt","w");
	double z = 0, thin = 10;
	for(z=0; z<=thin*(n-1);z++){
		fprintf(simple_out, "%10g %10g %10g\n",z/thin, akima_spline_eval(s,z/thin), cubic_spline_eval(cs,z/thin));
	}
	

	//doing a simple cos test
	int nn = 10;
	int i = 0;
	double xc[nn],yc[nn],dyc[nn];


	FILE* cospoint_out = fopen("cospoint_out.txt","w");
	for (i=0;i<nn;i++){
		xc[i] = 2*M_PI*i/nn;
		yc[i] = sin(xc[i]);
		dyc[i] = cos(xc[i]);
		fprintf(cospoint_out, "%10g %10g %10g\n",xc[i], yc[i], dyc[i]);
		}
	
//	akima_spline *ss = akima_spline_alloc(nn,xc,yc);
	cubic_spline *ccs = cubic_spline_alloc(nn,xc,yc,dyc);
	akima_spline *ss = akima_spline_alloc(nn,xc,yc);
	
	FILE* cos_out = fopen("cos_out.txt","w");
	double zz = 0.1; //I want to generate more points here then I tried for the last one
	for (double i =0; i<= 2*M_PI; i+=zz){
		double cubic_eval = cubic_spline_eval(ccs, i);
		fprintf(cos_out,"%10g %10g %10g\n",i,cubic_eval, akima_spline_eval(ss,i));
	}

	akima_spline_free(s);
	cubic_spline_free(cs);
	cubic_spline_free(ccs);

return 0;
}
