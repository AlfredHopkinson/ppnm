#include<stdio.h>
#include<math.h>

double (*f)(double); //says star f is a function and f is a pointer that takes a double then produces a function

void print_table(double f(double),double a, double b, double dx){
	for(double x=a;x<=b;x+=dx)
		printf("%10g %10g\n",x,f(x));
}


double square(double x){return x*x;}

int main(){
	double (*g)(double);
	g=&square;
	g=&cos;
	printf("g(2)=%g\n",g(2));
	print_table(g,0,M_PI,M_PI/8);
	print_table(sin,0,M_PI,M_PI/8);
return 0;
}
