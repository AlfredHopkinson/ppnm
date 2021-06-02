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


//activation function I have chosen is the fist one the gausian wavelet. I will change if the others seem better 
//double activation_func(double x){
//	return x*exp(-x*x);
//}

//void print_vector(const char* p, gsl_vector * A){
 //      	printf("%s\n",p);
//	for(int i=0;i<A->size;i++){
//			printf("%10g",gsl_vector_get(A,i));
//			}
//	printf("\n");
//}

//void linspace(gsl_vector * xs, double a, double b, int nn)  {
//
//	double dx = (b-a)/(nn-1);
//	for (int i = 0; i<nn; i++) {
//		gsl_vector_set(xs,i,a+i*dx);
//	}
//}
//
//double f(double x) {return sin(x);}












																
int main(){
	//allocate space for the network
	int n=3, nn = 20;
	double a = -2, b= 5;
	ann* network = ann_alloc(n, activation_func);
	
	for(int i=0;i<nn;i++){
		gsl_vector_set(network->params,3*i,a+(b-a)*i/(n-1));
		gsl_vector_set(network->params,3*i+1,1);
		gsl_vector_set(network->params,3*i+2,1);
		}

		//now I need to make the training and data   *****xs and ys remember******

		gsl_vector * xs = gsl_vector_alloc(nn);
		gsl_vector * ys = gsl_vector_alloc(nn);
		//had this going into out.txt but this got really messy so putting the output into num.txt

		FILE * num = fopen("num.txt","w");
		linspace(xs,a,b,nn);
		for(int i=0;i<nn;i++){
			gsl_vector_set(ys,i,f(gsl_vector_get(xs,i)));
			fprintf(num, "%10g %10g\n",gsl_vector_get(xs,i),gsl_vector_get(ys,i));
			}
		ann_train(network,xs,ys);
		print_vector("p = ",network->params);




	//remembering to free the data too
	
	ann_free(network);
	gsl_vector_free(ys);
	gsl_vector_free(xs);

return 0;
}





