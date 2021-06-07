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
#include"mini.h"
#include"ann.h"


void gradient(double F(gsl_vector * x), gsl_vector * x, gsl_vector * grad){
	double fx = F(x);
	int n = x->size;
	double DELTA = sqrt(DBL_EPSILON);
	for (int i = 0; i < n; i++) {
		double dx, ix = gsl_vector_get(x,i);
		if(fabs(ix)<DELTA) dx = DELTA;
		else dx = fabs(ix)*DELTA;
		gsl_vector_set(x,i,ix+dx);
		gsl_vector_set(grad,i,(F(x)-fx)/dx);
		gsl_vector_set(x,i,ix);
	}
}






int quasinewton(double F(gsl_vector*x), gsl_vector *x, double eps){
	int n = x->size;
	int iter =0;
	int goodstep=0,badstep=0;
	double DELTA = sqrt(DBL_EPSILON);

	gsl_matrix *H = gsl_matrix_alloc(n,n);
	gsl_vector *df = gsl_vector_alloc(n);
	gsl_vector *dx = gsl_vector_alloc(n);
	gsl_vector *prox = gsl_vector_alloc(n);
	gsl_vector *gs = gsl_vector_alloc(n);
	gsl_vector *y = gsl_vector_alloc(n);
	gsl_vector *u = gsl_vector_alloc(n);
	gsl_vector *step = gsl_vector_alloc(n);

	gsl_matrix_set_identity(H);
	gradient(F,x,df);
	double fiter;
	double fx = F(x);

	while(gsl_blas_dnrm2(df)>eps){
		iter++;
		gsl_blas_dgemv(CblasNoTrans,-1.0,H,df,0.0,dx);
		if(gsl_blas_dnrm2(dx) < DELTA*gsl_blas_dnrm2(x)) {fprintf(stderr,"|Dx|<DELTA*|x|\n"); break;} 
		if(gsl_blas_dnrm2(df) < eps) {fprintf(stderr,"|gradient|<eps\n"); break;} 
		double lambda = 1;

		while(1){
			gsl_vector_memcpy(step,x);
			gsl_vector_add(step,dx);
			double Hy;
			gsl_blas_ddot(dx,df,&Hy);
			fiter = F(step);
			if(fiter< fx+DELTA*Hy){goodstep++;break;}
			if(lambda < DELTA) { badstep++; gsl_matrix_set_identity(H); break; }
			lambda *=0.5;
			gsl_vector_scale(dx,0.5);
		}
		gradient(F,step,gs);
		gsl_vector_memcpy(y, gs);
		gsl_blas_daxpy(-1,df,y);
		gsl_vector_memcpy(u,dx);
		gsl_blas_dgemv(CblasNoTrans, -1,H,y,1.0,u);
		double T;
		gsl_blas_ddot(u,y,&T);
		if(fabs(T)>1e-10){
			gsl_blas_dger(1.0/T, u, u, H);
		}
		gsl_vector_memcpy(x,step);
		gsl_vector_memcpy(df,gs);
		fx=fiter;
	}
	gsl_matrix_free(H);
	gsl_vector_free(df);
	gsl_vector_free(dx);
	gsl_vector_free(prox);
	gsl_vector_free(gs);
	gsl_vector_free(y);
	gsl_vector_free(u);
	gsl_vector_free(step);

return iter;
}






//void gradient(void (*f)(gsl_vector * x, double * df), gsl_vector * x, gsl_vector * grad) {
//	double fx, faux; int n = x->size;
//	double DELTA = sqrt(DBL_EPSILON);
//	f(x,&fx);
//	for (int i = 0; i < n; i++) {
//		gsl_vector_set(x,i,gsl_vector_get(x,i)+DELTA);
//		f(x,&faux);
//		faux -= fx;
//		gsl_vector_set(grad,i,faux/DELTA);
//		gsl_vector_set(x,i,gsl_vector_get(x,i)-DELTA);
//	}
//}


//double backtrack(void (*f)(gsl_vector *x, double * fx), gsl_vector * x, gsl_vector * gf, gsl_vector * dx, gsl_matrix *H, int type){
//	assert(type == 1 || type ==0);
//
//	double l = 1, a = 0.0001, fx = f(x), deltaxdx;
//	gsl_vector *xandl = gsl_vector_alloc(x->size);
//	gsl_vector *ldx = gsl_vector_alloc((x->size));
//					
//	gsl_vector_memcpy(xandl,x);
//	gsl_vector_memcpy(ldx,dx);
//	gsl_vector_scale(ldx,l);
//	gsl_vector_add(xandl,ldx);
//	gsl_blas_ddot(dx,gf,&deltaxdx);
//	while( f(xandl)>= fx +a * l * deltaxdx){
//		l /= 2;
//		gsl_vector_memcpy(xandl,x);
//		gsl_vector_memcpy(ldx,dx);
//		gsl_vector_scale(ldx,l);
//		gsl_vector_add(xandl,ldx);
//
//		if(type==1){
//		if(l<0.000001){
//		gsl_matrix_set_identity(H);
//		break;}
//		}
//		else;
//	}
//	gsl_vector_free(xandl);
//	gsl_vector_free(ldx);
//
//	return l;
//}




//int quasinewton(double f(gsl_vector *x),gsl_vector *x,double eps){
//	int n = x->size;
//	int iter = 0;
//	
//	gsl_vector *df = gsl_vector_alloc(n);
//	gsl_vector *y = gsl_vector_alloc(n);
//	gsl_vector *s = gsl_vector_alloc(n);
//	gsl_matrix *H = gsl_matrix_alloc(n,n);
//	gsl_vector *Hy = gsl_vector_alloc(n);
//	double HyTs;
//
//	gsl_matrix_set_identity(H);
//	gradient(f,x,df);
//	while(gsl_blas_dnrm2(df)>eps){
//		iter++;
//		gsl_blas_dgemv(CblasNoTrans,-1.0,H,df,0.0,s);
//
//		double lambda = backtrack(f,x,df,s,H,1);
//
//		gsl_vector_scale(s,lambda);
//		gsl_vector_add(x,s);
//		gradient(f,x,y);
//		gsl_vector_sub(y,df);
//
//		gsl_blas_dgemv(CblasNoTrans,1.0,H,y,0.0,Hy);
//		gsl_blas_dgemv(CblasNoTrans,1.0,H,s,0.0,df);
//		gsl_blas_ddot(Hy,s,&HyTs);
//		gsl_vector_sub(s,Hy);
//		gsl_blas_dger(1.0/HyTs,s,df,H);
//
//		gradient(f,x,df);
//	}
//
//	gsl_vector_free(df);
//	gsl_vector_free(y);
//	gsl_vector_free(s);
//	gsl_matrix_free(H);
//	gsl_vector_free(Hy);
//return iter;
//}



//int quasinewton(double f(gsl_vector *x), void gradient(gsl_vector *x, gsl_vector *df),  gsl_vector *x,double eps){
//	int n = x->size;
//	int iter = 0;
//
//	gsl_vector *df = gsl_vector_alloc(n);
//	gsl_vector *y = gsl_vector_alloc(n);
//	gsl_vector *s = gsl_vector_alloc(n);
//	gsl_matrix *H = gsl_matrix_alloc(n,n);
//	gsl_vector *Hy = gsl_vector_alloc(n);
//	double HyTs;
//	gsl_matrix_set_identity(H);
//
//	gradient(x,df);
//
//	while( gsl_blas_dnrm2(df)>eps){
//		iter++;
//		gsl_blas_dgemv(CblasNoTrans,-1.0,H,df,0.0,s);
//
//		double lambda = backtrack(f,x,df,s,H,1);
//
//		gsl_vector_scale(s,lambda);
//		gsl_vector_add(x,s);
//		gradient(x,y);
//		gsl_vector_sub(y,df);
//		gsl_blas_dgemv(CblasNoTrans,1.0,H,y,0.0,Hy);
//		gsl_blas_dgemv(CblasNoTrans,1.0,H,s,0.0,df);
//		gsl_blas_ddot(Hy,s,&HyTs);
//		gsl_vector_sub(s,Hy);
//		gsl_blas_dger(1.0/HyTs,s,df,H);
//
//		gradient(x,df);
//	}
//
//	gsl_vector_free(df);
//	gsl_vector_free(y);
//	gsl_vector_free(s);
//	gsl_matrix_free(H);
//	gsl_vector_free(Hy);
//	return iter;
//}









