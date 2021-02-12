#include<stdio.h>
#include<limits.h>
#include<stdlib.h>
#include<assert.h>

//struct vector {int n, double a[]};
//struct vector my_vector;
//typedef struct ector vector

void set0(double x){x=0;}

void set0p(double *x){ (*x)=0; }

void print_array(int n, double a[]){
	for(int i=0;i<n;i++)printf("print_array: a[%d}=%g\n", i,i[a]);
	}

int main(){
	double y=1;
	set0(y);
	printf("y after set0 = %g\n", y);
	set0p(&y);
	printf("y after set0p = %g\n", y);
	
	int n=5;
	double v[5]; // define arrays like [] goes from 0to n-1
	for(int i=0;i<n;i++){
		v[i]=i;
		}	// i++ : i=i+1
	int i=0; while(i<n) {
		printf("v[%d]=%g\n",i,v[i]);
		i++;// this is a seperate i as the loop one only exists for the loop 
	}
	//if we put v[9999999]=1.0; it would be too large and push it over into a segmentation limit
	print_array(n,v);
	
	int N=7;
	double *a = malloc(N*sizeof(double)); //this asks a pointer to keep engough memeery allocated. its the same as the one above
	for(int i=0;i<N;i++) a[i]=i+100; //all of this memory is allocated manually and so has to be deleted
	print_array(N,a);

free(a);
return 0;
}
