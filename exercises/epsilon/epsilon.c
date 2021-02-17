#include<stdio.h>
#include<limits.h>
#include<stdlib.h>
#include<assert.h>
#include<float.h>

int main(){
	printf("Exercise 1 \n");
	int MAX = INT_MAX;
	int MIN = INT_MIN;

	printf("The max and min from limit.h are %i and %i\n", MAX, MIN);

	int iwhile = (MAX-10);
	int ifor = MAX-10;
	int ido = MAX-10;

	while(iwhile+1>iwhile){
		iwhile++;
		}
	for(int i=ifor;i+1>i; i++){
		ifor++;
	}
	do{ido++;} while( ido+1>ido);


	printf("My max with a while loop = %i\n",iwhile);
	printf("My max with a for loop = %i\n",ifor);
	printf("My max with a do while = %i\n",ido);

	printf("\n");
	
	int iwhilemin =1;
	int iformin =1;
	int idomin = 1;


	while(iwhilemin-1<1){
		iwhilemin--;
		}
	for(int i=0;i-1<1; i--){
		iformin--;
		}
	do{idomin--;} while(idomin-1<1);

	printf("My min with a while loop = %i\n",iwhilemin);
	printf("My min with a for loop = %i\n",iformin);
	printf("My min with a do while = %i\n",idomin);

	printf("\n");

	printf("Exercise iii\n");
	printf("epsilon machine while loops\n");

	double x=1;
	while(1+x!=1){
		x/=2;}
	x*=2;
	printf("double_epsilon = %g\n",x);
	printf("DBL_EPSILON = %g\n", DBL_EPSILON);

	long double ld=1;
	while(1+ld!=1){
		ld/=2;}
	ld*=2;
	printf("long_double_epsilon = %Lg\n",ld);
	printf("LDBL_EPSILON = %Lg\n",LDBL_EPSILON);

	float  fw = 1.0f;
	while(1+fw!=1){
		fw/=2;}
	fw*=2;
	printf("float_epsilon = %f\n",fw);
	printf("FLT_EPSILON = %f\n", FLT_EPSILON);

	printf("epsilon machine for loops\n");
	double x2=1;
	for(int i=0; 1+x2!=1;i++){x2/=2;} x2*=2;
	long double ld2=1;
	for(int i=0; 1+ld2!=1; i++) {ld2/=2;} ld2*=2;
	float f2=1;
	for(int i=0; 1+f2!=1; i++) {f2/=2;} f2*=2;

	printf("double_epsilon = %g\n",x2);
	printf("long_double_epsilon = %Lg\n",ld2);
	printf("float_epsilon = %f\n",f2);

	printf("do while loops for epsilon machine\n");
	double x3 = 1;
	do{x3/=2;} while(1+x3!=1); x3*=2;
	long double ld3 =1;
	do{ld3/=2;} while(1+ld3!=1); ld3*=2;
	float f3 =1;
	do{f3/=2;} while(1+f3!=1); f3*=2;
	printf("double_epsilon = %g\n",x3);
	printf("long_double_epsilon = %Lg\n",ld3);
	printf("float_epsilon = %f\n",f3);



	printf("\n");

	printf("Exercise 2\n");

	int ma = 7000000;
	int max = INT_MAX/3;

	float sum_up_float =0;
	float sum_down_float =0;

	for(int i=1;i<max;i++){sum_up_float += 1.0f/i;}
	printf("sum_up_float = %f\n", sum_up_float);

	for(int i=max;i>1;i--){sum_down_float += 1.0f/i;}
	printf("sum_down_float = %f\n",sum_down_float);

	printf("Exercise ii - They are different due to its float nature. it can only have a certain number of decimal points and so going big to small it adds same or smaller and this kleads to a difference.\n");


	printf("exercise iv \n");
	double sum_double_up = 0;
	double sum_double_down = 0;
	
	for(int i=1;i<max;i++){sum_double_up += 1.0f/i;}
	printf("sum_up_float = %f\n", sum_double_up);

	for(int i=max;i>1;i--){sum_double_down += 1.0f/i;}
	printf("sum_down_float = %f\n",sum_double_down);

	printf("Exercise 3\n");
	printf("This has been done in equal.c file");

	printf("here is an example a=30, b=30, tau=6, epsilon=6\n");
	int eq = equal(30,30,6,6);
	printf("this gives %i\n",eq);

	printf("here is an example a=25, b=35, tau=4, epsilon=4\n");
	int eq2 = equal(25,35,4,4);
	printf("this gives %i\n",eq2);
	




return 0;
}
