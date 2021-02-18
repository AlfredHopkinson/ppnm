#include<stdio.h>
#include<math.h>
#include <stdlib.h>

int main(int argc, char** argv){
	double x;
	int items;
	FILE* input=fopen("input.txt","r");
	FILE* out=fopen("out.file.txt","w");

	do{
		items=fscanf(input,"%lg",&x); // from stdin
		fprintf(out,"x=%g sin(x)=%g cos(x)=%g\n",x,sin(x),cos(x)); // to stdout
	}while(items!=EOF);


return 0;
}
