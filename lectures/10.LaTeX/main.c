#include<stdio.h>
#include<math.h>
double mygamma(double);
int main(){
	for(double x=1./8;x<=10;x+=1./8)
		printf("%g %g %g\n",x,mygamma(x),tgamma(x));
return 0;
}
