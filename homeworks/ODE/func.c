#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include<assert.h>

void print_vector(const char* p, gsl_vector * A){
	printf("%s\n",p);
	for(int i=0;i<A->size;i++){
		printf("%10g",gsl_vector_get(A,i));
	}
	printf("\n");
												
}
void print_matrix(const char* p, gsl_matrix * A){
	printf("%s\n",p);
	for(int i=0;i<A->size1;i++){
		for(int j=0;j<A->size2;j++){
			printf("%15g",gsl_matrix_get(A,i,j));
		}
		printf("\n");
	}
}

//make the RK stepper we go for 12 like in the example
void rkstep12(
	//compared to the exampel x=t, yx=yt, h=h, yh=yh, dy= err
	void f(double t,gsl_vector*y,gsl_vector*dydt) /* the f from dy/dt=f(t,y) */
	double t,              /* the current value of the variable */
	gsl_vector* yt,            /* the current value y(t) of the sought function */
	double h,              /* the step to be taken */
	gsl_vector* yh,             /* output: y(t+h) */
	gsl_vector* err){

	int c = yt->size;
	//dont need to alloc yt as its called above and alloced in main
	gsl_vector* k0 = gsl_vector_alloc(c);
	gsl_vector* k12 = gsl_vector_alloc(c);
	gsl_vector* ytt = gsl_vector_alloc(c);


	f(t,yt,k0); // this takes the step and evaluates 
	for(int i = 0; i<c; i++){
		gsl_vector_set(ytt,i, gsl_vector_get(yt,i)+gsl_vector_get(k0,i)*0.5*h);
	}
	f(t+0.5*h,ytt,k12);
	for(int i=0; i<c; i++){
		gsl_vector_set(yh,i, gsl_vector_get(yt,i)+gsl_vector_get(k12,i)*h);
	}
	for(int i=0;i<c;i++){
		gsl_vector_set(err,i, (gsl_vector_get(k0,i)-gsl_vector_get(k12,i))*0.5*h);
	}
	gsl_vector_free(k0);
	gsl_vector_free(k12);
	gsl_vector_free(ytt);

}

//now time to do the second driver part
//use same letters for the void like above
//follow the guide one in the chapter 
void driver(
	void f(double t,vector* y,vector* dydt), /* right-hand-side of dy/dt=f(t,y) */
	double a,                     /* the start-point a */
	gsl_vector* ya,                     /* y(a) */
	double b,                     /* the end-point of the integration */
	gsl_vector* yb,                     /* y(b) to be calculated */
	double h,                     /* initial step-size */
	double acc,                   /* absolute accuracy goal */
	double eps,
	double s,
	double s1,
	double err,
	double normyb,
	double tol
)
	while(t<b){
		if(t+h>b){ 
			h = b-t;//a check that it doesnt step too far
			rkstep12(f,t,ya,h,yb,err); //put it through the RK with the y(t=h) being now in yb
			
			//next part is checking the tolerance
			s=0;
			for(i=0; i<n; i++){ s+=gsl_vector_get(dy,i)*gsl_vector_get(dy,i);}
			err = sqrt(s);
			s1=0;
			for(i=0; i<n; i++){ s1+=gsl_vector_get(yb,i)*gsl_vector_get(yb,i);}
			normyb = sqrt(s1);
			tol = (normyb*eps+acc)*sqrt(h/(b-a));

			//this next part checks if the tolerance is greater then the error

			if(err<tol){
				t = t+h;
				gsl_vector_memcpy(ya,yb);
				for(int i =0; i<ya->size; i++){ printf("%g",gsl_vector_get(ya,i));
				}
				printf("%g",t);
			}
			h = pow(tol/err,o.25)*0.95; //the new stepsize
		}
		gsl_vector_free(yb);
		gsl_vector_free(err);
	}
