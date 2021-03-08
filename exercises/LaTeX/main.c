#include<stdio.h>
#include<math.h>

double exfun(double x){
	if(x<0)return 1/exfun(-x);
	if(x>1./8)return pow(exfun(x/2),2);
	return 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5*(1+x/6*(1+x/7*(1+x/8*(1+x/9*(1+x/10)))))))));
}

int main(){
	int x= -50;
	while(x<=50){
		printf("%10g %10g %10g\n",x/10.0,exfun(x/10.0),exp(x/10.0));
		x++;
	}
	
	return 0;
}
