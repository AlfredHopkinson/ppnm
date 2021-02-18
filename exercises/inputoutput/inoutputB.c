#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
	double x;
	int items;
	int i = 1;
	printf("Running B");
	do{
		items=fscanf(stdin,"%lg",&x);
		fprintf(stdout,"x=%g,sin(x)=%g,cos(x)=%g\n",x,sin(x),cos(x));
	}while(items!=EOF);
return 0;
}
