#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>
double t=0;
double fc (double x, void * params){
	
	double f = ((1/M_PI)*cos(x-t*sin(x)));
	return f;
}

double Bessel(double t){
	gsl_function F;
	F.function = &fc;
	int limit = 999;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (limit);
	double a=0, b=M_PI, acc = 1e-6, eps=1e-6, result, error;
	gsl_integration_qags(&F, a,b,acc,eps,limit,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}

int main(){
	for(double x=1;x<=20;x+=1.0/8){
		printf("%10g, %10g\n",x, Bessel(x));
		}
return 0;
}															
