#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include<assert.h>
#include"func.h"


//here make the SIR 
void sir(double t, gsl_vector* y, gsl_vector* dydt){
	//need pop size, time between contacts, recovery time
	double n = 5.83e6;
	double Tc = 3; //got this from sctoish covid 19 data
	double Tr = 14;
					
	//put in the pop dynamics
	gsl_vector_set(dydt,0, (-gsl_vector_get(y,1)*gsl_vector_get(y,0)/(n*Tc)));
	gsl_vector_set(dydt,1, (gsl_vector_get(y,1)*gsl_vector_get(y,0)/(n*Tc)-gsl_vector_get(y,1)/Tr));
	gsl_vector_set(dydt,2, (gsl_vector_get(y,1)/Tr));

}


int main(){


	double h =0.1;

	//initial parameters for Denamrk alloc 3 mnore matrices used the ya yb and error just as used in driver
	
	gsl_vector* ya = gsl_vector_alloc(3);
	gsl_vector* yb = gsl_vector_alloc(3);
	gsl_vector* error = gsl_vector_alloc(3);

	//initial conditions
	
	gsl_vector_set(ya,0,5.83e6);
	gsl_vector_set(ya,1,1e5); //had to raise this as I had it too high and it made it very fast
	gsl_vector_set(ya,2,0);

	//FILE * sirout = fopen("sir.out.txt","w");
	//driver for sir
	driver(sir,0,ya,60,yb,error,h,0.1,0.1);
	//fprintf(sirout,"%10d\n",sirsave);
	
	FILE * SIR_exp = fopen("SIR_Exp","w");
	printf("As can be seen in the graphs, increasing the Tc results in the infected peak becoming shorter and wider. It means the epidemic lasts longer because of this. We would expect this result with a longer time between infections.");



return 0;
}
