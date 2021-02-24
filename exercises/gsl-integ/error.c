#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

double fb (double x, void * params){
	double f = ((2/sqrt(M_PI)))*exp(-pow(x,2));
	return f;
}

double Berror(double x){
	gsl_function F;
	F.function = &fb;
	int limit = 999;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (limit);
	double a=0, b=x, acc = 1e-6, eps=1e-6, result, error;
	gsl_integration_qags(&F, a,b,acc,eps,limit,w,&result,&error);
	gsl_integration_workspace_free(w);
	return result;
}




int main(){
	for(double x=-3;x<=3;x+=1.0/8){
		printf("%10g, %10g\n",x, Berror(x));
		}
return 0;
}															
