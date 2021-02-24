#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

 
double fa (double x, void * params){
	double f = (log(x)/sqrt(x));
	return f;
}

double Aint(){
	gsl_function F;
	F.function = &fa;
	int limit = 999;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (limit);
	double a=0, b=1, acc = 1e-6, eps=1e-6, result, error;
	gsl_integration_qags(&F, a,b,acc,eps,limit,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}

int main(){
	printf("When using gsl qags the value of the intergration is %g\n",Aint());
	
return 0;
}															
