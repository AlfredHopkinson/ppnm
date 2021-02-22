#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>
#include<gsl/gsl_sf_gamma.h>

double Erf(double);
double Gamma(double);


int main(){
	FILE* data_file = fopen("data.txt","w");
	double xmin =-2, xmax=2;
	for(double x=xmin; x<=xmax; x+=1.0/8){
		fprintf(data_file, "%10g %10g %10g %10g\n",x,erf(x),gsl_sf_erf(x),Erf(x));
	}
	fclose(data_file);

	FILE* gamma_file = fopen("gamma.txt","w");
	double gamin = -2, gamax =5;
	for(double x=gamin; x<=gamax; x+=1.0/9){
		fprintf(gamma_file," %10g %10g %10g\n",x,tgamma(x), Gamma(x));
	}
	fclose(gamma_file);


return 0;
}
