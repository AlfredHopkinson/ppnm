#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>
#include<gsl/gsl_sf_gamma.h>
#define M_PI 3.1415926535897

double Gamma(double x){
	///single precision gamma function (Gergo Nemes, from Wikipedia)
	if(x<0)return M_PI/sin(M_PI*x)/Gamma(1-x);
	if(x<9)return Gamma(x+1)/x;
	double lnGamma=x*log(x+1/(12*x-1/x/10))-x+log(2*M_PI/x)/2;
	return exp(lnGamma);
	}

void make_gamma_data(FILE* file)
{
	double xmin =-2, xmax=5;
	for(double x=xmin; x<=xmax; x+=1.0/8){
		printf("%10g %10g %10g %10g\n",x, tgamma(x), Gamma(x), log(tgamma(x)));
	}
}




int main(){
	FILE* gammaf = fopen("gamma.txt","w");
	make_gamma_data(gammaf);


return 0;
		
}
