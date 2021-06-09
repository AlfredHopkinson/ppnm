#include<stdio.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
#include<assert.h>
#include <gsl/gsl_integration.h>
#define RND ((double)rand()/RAND_MAX)
#include<float.h>
#include"ann.h"
#include"mini.h"


int n =3;
ann* network;
int x = 20;
gsl_vector*xs;
gsl_vector*ys;

//activation function I have chosen is the fist one the gausian wavelet. I will change if the others seem better 
double activation_func(double x){
	return x*exp(-x*x);
}

void print_vector(const char* p, gsl_vector * A){
      	printf("%s\n",p);
	for(int i=0;i<A->size;i++){
			printf("%10g",gsl_vector_get(A,i));
			}
	printf("\n");
}

void linspace(gsl_vector * xs, double a, double b, int nn)  {

	double dx = (b-a)/(nn-1);
	for (int i = 0; i<nn; i++) {
		gsl_vector_set(xs,i,a+i*dx);
	}
}
//function to fit 
double f(double x){return sin(x);} 



double cf (gsl_vector * c){
	gsl_vector_memcpy(network->params,c);
	double count =0;
	for (int i=0; i<xs->size; i++){
		double ix = gsl_vector_get(xs,i);
		double iy = gsl_vector_get(ys,i);
		double fi = ann_response(network,ix);
		count += fabs(fi-iy);
	}
	return count/xs->size;
}



//for future me I had a lot of difficulty here and then debuging was vbery hard. the spaces serve no purpose it was where I had print functions to make it easy to see the fault. if using for another thing then just ignore these



																
int main(){
	//allocate space for the network
	int n=3, nn = 50;
	double a = -3.14, b= 3.14;
//	printf("PART 1\n");
	//ann* network = ann_alloc(n, activation_func);
	network = ann_alloc(n,activation_func);

//	printf("PART 2\n");
//	assign first params
//	for(int i=0;i<n;i++){
//		gsl_vector_set(network->params,3*i,a+(b-a)*i/(n-1));
//		gsl_vector_set(network->params,3*i+1,1);
//		gsl_vector_set(network->params,3*i+2,1);
//		
//		}
		

		//now I need to make the training and data   *****xs and ys remember******
//		printf("PART 4\n");
		//gsl_vector * xs = gsl_vector_alloc(nn);
		//gsl_vector * ys = gsl_vector_alloc(nn);
		xs = gsl_vector_alloc(nn);
		ys = gsl_vector_alloc(nn);

		//generating x and y numbers
	for(int i=0; i<nn; i++){
		double x = a+(b-a)*i/(nn-1);
		double fy = activation_func(x);
		gsl_vector_set(xs,i,x);
		gsl_vector_set(ys,i,fy);
	}

	for(int i=0;i<n;i++){
		gsl_vector_set(network->params,3*i,a+(b-a)*i/(n-1));
		gsl_vector_set(network->params,3*i+1,1);
		gsl_vector_set(network->params,3*i+2,1);
	}



		//had this going into out.txt but this got really messy so putting the output into num.txt
//		printf("PART 5\n");
		FILE * num = fopen("num.txt","w");
//		linspace(xs,a,b,nn);
//		printf("PART 6\n");
		for(int i=0;i<nn;i++){
			gsl_vector_set(ys,i,f(gsl_vector_get(xs,i)));
			fprintf(num, "%10g %10g\n",gsl_vector_get(xs,i),gsl_vector_get(ys,i));
			}
//		printf("PART 7\n");
		ann_train(network,xs,ys);
//		printf("PART 8\n");
//		print_vector("p = ",network->params);

		for(int i=0; i<nn; i++) {
			double x = gsl_vector_get(xs, i);
			double f = gsl_vector_get(ys, i);
			double val = ann_response(network,x);
			printf("%10g %10g %10g\n", x, f,val);
		}
											




	//remembering to free the data too
	
	ann_free(network);
	gsl_vector_free(ys);
	gsl_vector_free(xs);

return 0;
}





