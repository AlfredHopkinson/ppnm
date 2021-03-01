#include<stdio.h>
#include<pthread.h>
#include<math.h>
#include<stdlib.h>

#define RND (double)rand_r()/RAND_r_MAX



void* bar(void* arg){
	double* incircle = (double*)arg;
	double x,y;
	unsigned int seedx = 1;
	unsigned int seedy = 2;

	for(int i=0;i<1e8;i++){
		x = rand_r(&seedx);
		y = rand_r(&seedy);
		x/=RAND_MAX;
		y/=RAND_MAX;
		if (pow(x,2)+pow(y,2)<=1){
			*incircle=*incircle+1;
		}
	}
				
return NULL;
}

int main(){
	

	double x=0, y=0, z=0;
	pthread_t threadx, thready;
	pthread_attr_t* attributes = NULL;
	pthread_create( &threadx, attributes, bar, (void*)&x);
	pthread_create( &thready, attributes, bar, (void*)&y);
							
	bar((void*)&z); //meanwhile in the main thread

	void* returnvalue = NULL;
	pthread_join(threadx,returnvalue);
	pthread_join(thready,returnvalue);

	double pi = 4*((x+y+z)/3e8);
	
	
	printf("Value of Pi is %g\n",pi);
	
return 0;
}




//define the square by an array
//put in by a random functions
//use bar to get number to throw as a paramter and gives number of points in a circle back 
//call it several times then do the calc
//use x + y both squared top get the circle and if its greater than 1 its out
//
//
//create funt to make points and checks if inside then gives total 
//run this in multiple threads
