#include<math.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>

double lin_vec_interp(gsl_vector x, gsl_vector y, double z);
double lin_vec_integ(gsl_vector x, gsl_vector y, double z);
double z;

int main(){
	
	int n = 15;
	gsl_vector *x = gsl_vector_alloc(n+1);
	gsl_vector *y = gsl_vector_alloc(n+1);
	for(int i=0;i<=n;i++)
		{
		gsl_vector_set (x, i, i);
		gsl_vector_set (y, i, 1);
		}

	
	FILE * linfile=fopen("lin.out.txt","w");
	for(z=1.0/8;z<15;z+=1.0/8){
	double interp = linterp(x,y,z);

	fprintf(linfile,"%g %g\n",z,interp);}

fclose(linfile);

gsl_vector_free(x);
gsl_vector_free(y);
return 0;
}


