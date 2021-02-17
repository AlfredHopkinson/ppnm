#include <stdio.h>
#include<math.h>

void f(int* i){*i=0;}

int main(void){
	double x=1.23;
	double y = *(&x);
	printf("x = %g\n",x);
	printf("y = %g\n",y);

	int i=1; f(&i); printf("i=%i\n",i);

return 0;
}
