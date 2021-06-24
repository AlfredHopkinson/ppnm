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

//make the gradient it says
//void gradient (void(*f)(gsl_vector *x, double fx), gsl_vector *gradianfx, gsl_vector *dx){
//	double l = 1, in;
//	double a = 0.000001
//	double n = 0, maxn=1000;
//	double fx, fxs;
//	gsl_vector * xl = gsl_vector_alloc(x->size);
//	while (n<maxn){
//		n++;
//		gsl_vector_memcpy(xl,x);
//		gsl_vector_add(xl,dx);
//		
//		(x,&fx);
//		gsl_blas_ddot(dx,gradianfx,&in);
//		if (l < DELTA) break;
//		if (fxs < fx + a*in) break;
//		l /= 2.0;
//		gsl_vector_scale(dx,l);
//	}
//	gsl_vector_free(xl);
//}


//we need to make the quasinewton func and have the implimentatiopn template to work off.
//double DELTA = sqrt(DBL_EPSILON);
//void qnewton(void (*f)(gsl_vector * x, gsl_vector * fx),gsl_vector * x, double eps){
//	int n = x->size, numsteps = 0;
//	
//	gsl_vector * gfx = gsl_vector_alloc(n);
//	gsl_vector * e = gsl_vector_alloc(n);
//
//	while (numsteps<1000){
//		if(gsl_blas_dnrm2(gfx)<eps*eps) break;
//		gsl_blas_dgemv(CblasNoTrans,-1,B,gfx,0,e);
//		gradient(f,x,gfx,e);

//add in the rosenbrock stuff need it, its gradient and hessian

double rosenbrock(gsl_vector *x) 
{
	double xa = gsl_vector_get(x,0);
	double ya = gsl_vector_get(x,1);

	return pow(1-xa,2)+100*pow(ya-pow(xa,2),2);
}

void rosenbrockgrad (gsl_vector * R , gsl_vector * fR){	
	double gx,gy;
	double xa = gsl_vector_get(R,0);
	double ya = gsl_vector_get(R,1);
	gx = -2*(1-xa)-400*xa*(ya-(xa*xa));
	gy = 200*(ya-(xa*xa));
	gsl_vector_set(fR,0,gx);
	gsl_vector_set(fR,1,gy);
}

void rosenbrockhess(gsl_vector *x, gsl_matrix *H) 
{
	double xa = gsl_vector_get(x,0);
	double ya = gsl_vector_get(x,1);

	gsl_matrix_set(H,0,0,2-400*ya+1200*pow(xa,2));
	gsl_matrix_set(H,1,1,200);
	gsl_matrix_set(H,0,1,-400*xa);
	gsl_matrix_set(H,1,0,-400*xa);
}

//him stuff here

double himmel(gsl_vector *x)
{
	double xa = gsl_vector_get(x,0);
	double ya = gsl_vector_get(x,1);

	return pow(pow(xa,2)+ya-11,2)+pow(xa+pow(ya,2)-7,2);
}

void himmelgrad(gsl_vector *h ,gsl_vector *fh)
{
	double gx,gy;
	double xa = gsl_vector_get(h,0);
	double ya = gsl_vector_get(h,1);
	gx = 4*xa*(pow(xa,2)+ya-11)+2*(xa+pow(ya,2)-7);
	gy = 2*(pow(xa,2)+ya-11)+4*ya*(xa+pow(ya,2)-7);
	gsl_vector_set(fh,0,gx);
	gsl_vector_set(fh,1,gy);
}

void himmelhess(gsl_vector *x, gsl_matrix *H)
{
	double xa = gsl_vector_get(x,0);
	double ya = gsl_vector_get(x,1);
	double H12 = 4*(xa+ya);
	gsl_matrix_set(H,0,0,4*ya+12*pow(xa,2)-2);
	gsl_matrix_set(H,0,1,H12);//H12 and 21 are the same
	gsl_matrix_set(H,1,0,H12);
	gsl_matrix_set(H,1,1,4*xa+12*pow(ya,2)-28);
}






//Make the gradient

void gradient(void (*f)(gsl_vector * x, double * df), gsl_vector * x, gsl_vector * grad) {
	double fx, faux; int n = x->size;
	double delt = sqrt(DBL_EPSILON);
	f(x,&fx);
	for (int m = 0; m < n; m++) {
		gsl_vector_set(x,m,gsl_vector_get(x,m)+delt);
		f(x,&faux); faux -= fx;
		gsl_vector_set(grad,m,faux/delt); gsl_vector_set(x,m,gsl_vector_get(x,m)-delt);
	}
}


//backtrack here

double backtrack(double f(gsl_vector *x), gsl_vector * x, gsl_vector * gf, gsl_vector * dx, gsl_matrix *H, int type){
	//here put in the broyden update
	assert(type == 1 || type ==0);

	double l = 1, a = 0.0001, fx = f(x), deltaxdx;
	gsl_vector *xandl = gsl_vector_alloc(x->size);
	gsl_vector *ldx = gsl_vector_alloc((x->size));
	
	gsl_vector_memcpy(xandl,x);
	gsl_vector_memcpy(ldx,dx);
	gsl_vector_scale(ldx,l);
	gsl_vector_add(xandl,ldx);

	gsl_blas_ddot(dx,gf,&deltaxdx);

	while( f(xandl)>= fx +a * l * deltaxdx){
		l /= 2;
		gsl_vector_memcpy(xandl,x);
		gsl_vector_memcpy(ldx,dx);
		gsl_vector_scale(ldx,l);
		gsl_vector_add(xandl,ldx);

		if(type==1){
		if(l<0.000001){
		gsl_matrix_set_identity(H);
		break;}
		}
		else;
	}
	gsl_vector_free(xandl);
	gsl_vector_free(ldx);

	return l;
}



	
int quasinewton(double f(gsl_vector *x),void gradient(gsl_vector *x, gsl_vector *df),gsl_vector *x,double eps){
	int n = x->size;
	int iter = 0;

	gsl_vector *df = gsl_vector_alloc(n);
	gsl_vector *y = gsl_vector_alloc(n);
	gsl_vector *s = gsl_vector_alloc(n);
	gsl_matrix *H = gsl_matrix_alloc(n,n);
	gsl_vector *Hy = gsl_vector_alloc(n);
	double HyTs;
	gsl_matrix_set_identity(H);

	gradient(x,df);

	while( gsl_blas_dnrm2(df)>eps){
		iter++;
		gsl_blas_dgemv(CblasNoTrans,-1.0,H,df,0.0,s);

		double lambda = backtrack(f,x,df,s,H,1);

		gsl_vector_scale(s,lambda);
		gsl_vector_add(x,s);
		gradient(x,y);
		gsl_vector_sub(y,df);

		gsl_blas_dgemv(CblasNoTrans,1.0,H,y,0.0,Hy);
		gsl_blas_dgemv(CblasNoTrans,1.0,H,s,0.0,df);
		gsl_blas_ddot(Hy,s,&HyTs);
		gsl_vector_sub(s,Hy);
		gsl_blas_dger(1.0/HyTs,s,df,H);

		gradient(x,df);
	}

	gsl_vector_free(df);
	gsl_vector_free(y);
	gsl_vector_free(s);
	gsl_matrix_free(H);
	gsl_vector_free(Hy);
	return iter;
}








int main(){
	int n = 2, iter;
	gsl_vector *x = gsl_vector_alloc(n);
	gsl_vector *df = gsl_vector_alloc(n);
	double eps = 1e-6;

	//rosen part A
	printf("PART A \n\n\n");

	printf("Rosenbrock\n\n");

	gsl_vector_set(x,0,5);
	gsl_vector_set(x,1,10);

	iter = quasinewton(&rosenbrock,&rosenbrockgrad,x,eps);

	rosenbrockgrad(x,df);

	printf("rosen min is x = %10g y = %10g\n",gsl_vector_get(x,0), gsl_vector_get(x,1));
	printf("num iterations = %10d\n",iter);

	//himmel part a
	printf("Himmel\n\n");

	gsl_vector_set(x,0,5);
	gsl_vector_set(x,1,10);

	iter = quasinewton(&himmel, &himmelgrad,x,eps);

	himmelgrad(x,df);

	printf("himmel min is x = %10g y = %10g\n",gsl_vector_get(x,0), gsl_vector_get(x,1));
	printf("num iterations = %10d\n",iter);
	


return 0;
}









