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

//as said in the lecture here i can use structs to keep it in so i dont calculate each time
// this was in the booklet and so i used it here for the quadratics


typedef struct {int n; double *x, *y, *b, *c;} qspline;

qspline * qspline_alloc(int n, double *x, double *y){ //this is returning a pointer to the structure
	qspline *s = (qspline*)malloc(sizeof(qspline));
	s->c = (double*) malloc((n-1)*sizeof(double));
	s->b = (double*) malloc((n-1)*sizeof(double));
	s->x = (double*) malloc(n*sizeof(double));
	s->y = (double*) malloc(n*sizeof(double));
	s->n = n;
	int i;
	for(i=0;i<n;i++){s->x[i]=x[i];s->y[i]=y[i];}
	double p[n-1], h[n-1];
	for(i=0;i<n-1;i++){h[i]=x[i+1]-x[i];p[i]=(y[i+1]-y[i])/h[i];}
	s->c[0]=0;
	for (i=0;i<n-2;i++) {s->c[i+1] = (p[i+1]-p[i]-s->c[i]*h[i])/h[i+1];}
	s->c[n-2]/=2;
	for (i=n-3;i>=0;i--) {s->c[i] = (p[i+1]-p[i]- s->c[i+1]*h[i+1])/h[i];}
	for(i=0;i<n-1;i++) {s->b[i]=p[i]-s->c[i]*h[i];}
	return s;
}
	
double qspline_eval(qspline *s, double z){
	assert (z>=s->x[0] && z<=s->x[s->n-1]);
	int i=0, j=s->n-1;
	while(j-i>1){int m=(i+j)/2; if (z>s->x[m]) i=m; else j=m;}
	double h=z-s->x[i];
	return s->y[i]+h*(s->b[i]+h*s->c[i]);}


//do the derivitive
double qspline_derivitive(qspline *s, double z){
	assert (z>=s->x[0] && z<=s->x[s->n-1]); //allways keep the assert rules in and take i as the binsearch is easiest  remember to swtich out all others where you did it from the chapter going forward
	int i = binsearch(s->n,s->x,z);
	double dx = z - s->x[i];
	return s->b[i] + 2*s->c[i]*dx;
}

//do the integral

double qspline_integ(qspline *s, double z){
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

//had to put this at the end as it didnt work for ages up at the top as this frees the allocated memory
void qspline_free(qspline *s){
	free(s->x); free(s->y); free(s->b); free(s->c); free(s);}


//do the main
int main(){
	//need the allocation of the points 
	int n=10;
	double x[n], y[n];
	FILE* mypoints = fopen("plots.out.txt","w");

	for(int i=0;i<n;i++){
		x[i]=i;
		y[i]=pow(i,2);
		fprintf(mypoints,"%10g %10g\n",x[i],y[i]);
	}

	double thin = 10;
	int z = 0; //taken and modified from linear
	qspline * s = qspline_alloc(n,x,y);
	
	//insead of many many different files put all into one file 

	FILE* quadout=fopen("quad.out.txt","w");
	for(z=0; z<=thin*(n-1);z++){
		fprintf(quadout,"%10g %10g %10g %10g\n",z/thin, qspline_eval(s,z/thin), qspline_derivitive(s,z/thin), qspline_integ(s,z/thin));
	}
return 0;		
}






	
