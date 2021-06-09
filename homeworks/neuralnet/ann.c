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
//using the same layoput for this file and the .h file as suggested in the question



ann* ann_alloc(int n, double(*f)(double)){

	ann * network = malloc(sizeof(ann));
	network->n = n;
	network->f = f;
	network->params = gsl_vector_alloc(3*n);
	return network;

}


void ann_free(ann* network){
	//free the netwoork
	gsl_vector_free(network->params);
	free (network);


}

double ann_response (ann* network, double x){

	double response = 0;
	for(int i = 0; i<network->n; i++){//set the paramters
		double a = gsl_vector_get(network->params,0*network->n+i);
		double b = gsl_vector_get(network->params,1*network->n+i);
		double w = gsl_vector_get(network->params,2*network->n+i);
		response += network -> f((x+a)/b)*w;
	}
	return response;

}

double cf(gsl_vector* q);

void ann_train (ann* network, gsl_vector* xs, gsl_vector* ys){

//	void delta(gsl_vector * p)
//	{
//		gsl_vector_memcpy(network->params,p);
//		double response = 0;
//		for(int i=0; i<xs->size; i++){
//			double x = gsl_vector_get(xs,i);
//			double y = ann_response(network,x);
//			double f = gsl_vector_get(ys,i);
//			response += fabs(y-f);
//		}
//		return response/xs -> size;
//	}


	
	gsl_vector *q = gsl_vector_alloc(network->params->size);
	gsl_vector_memcpy(q,network->params);
	quasinewton(cf,q,1e-3);
	gsl_vector_memcpy(network->params,q);
	
	gsl_vector_free(q);
}


	
//	gsl_vector * q = gsl_vector_alloc(3*network->n);
//	gsl_vector_memcpy(q, network->params);
//	int step = quasinewton(delta,q,1e-3);
//	gsl_vector_memcpy(network->params,q);
//	printf("The number of steps taken to train is = %i\n",step);
//	gsl_vector_free(q);
	//do the cost function
	



	//double s = 0.0001;
	//int step = quasinewton(delta,gsl_vector * q,s);
	//gsl_vector_memcpy(network->params, q);
	//gsl_vector_free(q);
	

//}



